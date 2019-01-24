// parse barcode sequence from raw sequences for sci library and QC
#include "utils.h"
#include "fastq.h"
#include "thread_pool.h"
#include "number.h"
#include <zlib.h>
#include <htslib/kstring.h>
#include <ctype.h>

#define PLT_UNKNOWN -1
#define PLT_SCI      0
#define PLT_DROP     1
#define PLT_RENAME   2
#define PLT_STREAM   3 // do nothing, just streaming the fastqs
#define PLT_PLATE    4

static char *program_name =  "sc_atac_parse_barcode";

static char *version = "0.0.0.9002";

static int get_version()
{
    puts(version);
    return 1;
}

int usage()
{
    fprintf(stderr, "* %s read_1.fq.gz read_2.fq.gz\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -1 out_r1.fq                         R1 output.\n");
    fprintf(stderr, "  -2 out_r2.fq                         R2 output.\n");
    fprintf(stderr, "  -report split_report.txt             Process report, in markdown format.\n");
    fprintf(stderr, "  -cb cell_barcode.txt                 Cell barcode sequence and quality score.\n");
    fprintf(stderr, "  -platform [sci|drop|plate|stream]    Select platform to process.\n");
    fprintf(stderr, "  -lane                                Set the lane number.\n");
    fprintf(stderr, "  -t [5]                               Threads.\n");    
    fprintf(stderr, "  -r [100000]                          Records per chunk.\n");
    fprintf(stderr, "  -p                                   Read 1 and read 2 in same file.\n");
    fprintf(stderr, "  -name CELL_BARCODE                   Manually set the cell barcode for the fastq(s), skip parse from the sequence.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Notes: default output smart read pairs in stdout.\n");
    return 1;
}

struct args {
    const char *r1_fname;
    const char *r2_fname;
    const char *out1_fname;
    const char *out2_fname;
    const char *report_fname;
    const char *barcode_fname;

    const char *name;

    const char *lane;
    gzFile r1_fp;
    gzFile r2_fp;
    FILE *out1_fp;
    FILE *out2_fp;
    FILE *report_fp;
    FILE *barcode_fp;
    
    const char *platform;
    int platform_id;
    
    int rm_qc_failed; // qc failed reads will not be exported
    int rm_unknown_bc; // unknown barcodes will not be exported

    int n_thread;
    int chunk_size;

    // flag set for input, program will auto identify output mode
    int smart_pair;
    //int is_pe;
    
    struct fastq_handler *fastq;
    struct qc_report report;
} args = {
    .r1_fname = NULL,
    .r2_fname = NULL,
    .out1_fname = NULL,
    .out2_fname = NULL,
    .report_fname = NULL,
    .barcode_fname = NULL,

    .name = NULL,
    .lane = NULL,
    .r1_fp = NULL,
    .r2_fp = NULL,
    .out1_fp = NULL,
    .out2_fp = NULL,
    .report_fp = NULL,
    .barcode_fp = NULL,
    
    .platform = NULL,
    .platform_id = PLT_UNKNOWN,
    
    .rm_qc_failed = 1,
    .rm_unknown_bc = 1,

    .n_thread = 5,
    .chunk_size = 100000,

    .smart_pair = 0,
    //.is_pe = -1,

    .fastq = NULL,
    .report = {0,0,0},
};

static int check_baseN(char *s, int l)
{
    int i;
    for (i = 0; i < l; ++i) {
        if (s[i] == 'N') return 1;
    }
    return 0;
}
// add cell barcode name directly, no QC, no parsing
void *manually_rename(void *_data, int idx)
{
    struct bseq_pool *p = (struct bseq_pool *)_data;
    struct args *opts = (struct args*)p->opts;
    assert(opts->name);
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        assert(b && b->n0);
        // it is silly to export same cell barcode name for each read, but we just make this to compare our pipeline
        kstring_t cb = {0,0,0};
        kputs(opts->name, &cb); 
        b->data = (void*)cb.s;
        
        kstring_t str = {0,0,0};
        kputs(b->n0, &str);
        kputs("|||CR|||", &str);
        kputs(opts->name, &str);        
        free(b->n0);
        b->n0 = str.s;
        b->flag = FQ_PASS;
    }
    return (void*)p;
}
void *stream_run(void *_data, int idx)
{
    struct bseq_pool *p = (struct bseq_pool *)_data;
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        b->flag = FQ_PASS;

        kstring_t cb = {0,0,0};
        char *s = b->n0;
        char *e = s+(strlen(s)-1);
        while(s && *s != '|' && s != e) s++;
        if (s) { // *s == '|' definitely
            char *p = s;
            int j = 0;
            while(p && p != e) {
                p++;
                j++;
            }
            if (j>8) {
                if (s[1] == '|' && s[2] == '|' && s[3] == 'C' && s[4] == 'R' && s[5] == '|' && s[6] == '|' && s[7] == '|') {
                    // CR found
                    s = s + 8;
                    p = s;
                    int k = 0;
                    while (p != e && *p != '|' && !isspace(*p)) { p++; k++; }
                    kputsn(s, k, &cb);
                    b->data = (void*)cb.s;
                }
            }
        }
    }
    return _data;
}
void *plate_run(void *_data, int idx)
{
    struct bseq_pool *p = (struct bseq_pool *)_data;
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        assert (b->l0 == 50 && b->l1 == 58); // predefined format
        // TODO: check barcode list
        kstring_t cb = {0,0,0};
        kputsn(b->s1+50, 8, &cb);

        if ( check_baseN(cb.s, cb.l)) b->flag = FQ_BC_FAIL;

        kputc('\t', &cb);
        kputsn(b->q1+50, 8, &cb);
        kputs("", &cb);
        b->data = (void*)cb.s;
           
        // read structure may be updated in the futher
        kstring_t str = {0,0,0};
        kputs(b->n0, &str);
        kputs("|||CR|||", &str);
        kputsn(b->s1+50, 8, &str);        
        kputs("|||CY|||", &str);
        kputsn(b->q1+50, 8, &str);
        kputs("", &str);

        b->s1[50]='\0';
        b->q1[50]='\0';
        b->l1 = 50;

        free(b->n0);
        b->n0 = str.s;
    }
    return (void*)p;
}

void *sci_run(void *_data, int idx)
{
/*
  The following format is supplied by Wu Liang.
  
  $b1=substr($r1_2,0,10);
  $b2=substr($r1_2,31,10);
  $s1=substr($r1_2,60,34);
  $b3=substr($r2_2,0,10);
  $b4=substr($r2_2,37,10);
  $s2=substr($r2_2,66,34);
  $q1=substr($r1_4,60,34);
  $q2=substr($r2_4,66,34);
  print OR1 "\@$b1$b2$b3$b4:$num\n$s1\n+\n$q1\n";
  print OR2 "\@$b1$b2$b3$b4:$num\n$s2\n+\n$q2\n";
*/
    
    struct bseq_pool *p = (struct bseq_pool *)_data;
    struct args *opts = (struct args*)p->opts;
    
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        assert (b->l0 == 100 && b->l1 == 100); // predefined format
        // TODO: check barcode list
        kstring_t cb = {0,0,0};
        kputsn(b->s0, 10, &cb);
        kputsn(b->s0+31, 10, &cb);
        kputsn(b->s1,10, &cb);
        kputsn(b->s1+37, 10, &cb);

        if ( check_baseN(cb.s, cb.l)) b->flag = FQ_BC_FAIL;
             
        kputc('\t', &cb);
        kputsn(b->q0, 10, &cb);
        kputsn(b->q0+31, 10, &cb);
        kputsn(b->q1, 10, &cb);
        kputsn(b->q1+37, 10, &cb);
        kputs("", &cb);
        b->data = (void*)cb.s;
           
        // read structure may be updated in the futher
        kstring_t str = {0,0,0};
        kputs(b->n0, &str);
        kputs("|||CR|||", &str);
        kputsn(b->s0, 10, &str);
        kputsn(b->s0+31,10, &str);
        kputsn(b->s1,10, &str);
        kputsn(b->s1+37,10, &str);
        if (opts->lane) {
            kputc('-', &str);
            kputs(opts->lane, &str);
        }
        kputs("|||CY|||", &str);
        kputsn(b->q0, 10, &str);
        kputsn(b->q0+31,10, &str);
        kputsn(b->q1,10, &str);
        kputsn(b->q1+37,10, &str);

        memmove(b->s0, b->s0+60, sizeof(char)*34);
        memmove(b->q0, b->q0+60, sizeof(char)*34);
        b->s0[34]='\0';
        b->q0[34]='\0';
        b->l0 = 34;
       
        memmove(b->s1, b->s1+66, sizeof(char)*34);
        memmove(b->q1, b->q1+66, sizeof(char)*34);
        b->s1[34]='\0';
        b->q1[34]='\0';
        b->l1 = 34;

        free(b->n0);
        b->n0 = str.s;
    }
    return (void*)p;
}

void *drop_run(void *_data, int idx)
{
    struct bseq_pool *p = (struct bseq_pool *)_data;
    int i;
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        if (b->l1 == 0) continue; // if no read 2, skip
        if (b->l1 == 88) { // predefined format
            // TODO: check barcode list
            kstring_t cb = {0,0,0};            
            kputsn(b->s1+50,10, &cb);
            kputsn(b->s1+78,10, &cb);

            if ( check_baseN(cb.s, cb.l)) b->flag = FQ_BC_FAIL;
            
            kputc('\t', &cb);
            kputsn(b->q1+50,10, &cb);
            kputsn(b->q1+78,10, &cb);
            kputs("", &cb);
            b->data = (void*)cb.s;
            
            // read structure may be updated in the futher
            kstring_t str = {0,0,0};
            kputs(b->n0, &str);
            kputs("|||CR|||", &str);
            kputsn(b->s1+50,10, &str);
            kputsn(b->s1+78,10, &str);
            kputs("|||CY|||", &str);
            kputsn(b->q1+50,10, &str);
            kputsn(b->q1+78,10, &str);
            b->s1[50]='\0';
            b->q1[50]='\0';
            b->l1 = 50;
            //
            free(b->n0);
            b->n0 = str.s;
        }
        else error("Unsupport format. %d", b->l1);
    }
    return (void*)p;
}
void write_out(void *_data)
{
    struct bseq_pool *p = (struct bseq_pool*)_data;
    struct args *opts = (struct args*)p->opts;

    FILE *fp1 = opts->out1_fp == NULL ? stdout : opts->out1_fp;
    FILE *fp2 = opts->out2_fp == NULL ? fp1 : opts->out2_fp;
    int i;
    // because the output queue is order, we do not consider the thread-safe of summary report
    for (i = 0; i < p->n; ++i) {
        struct bseq *b = &p->s[i];
        opts->report.all_fragments++;
        
        if (b->flag == FQ_QC_FAIL) {
            opts->report.qc_failed++;
            continue;
        }
        else if (b->flag == FQ_BC_FAIL) {
            opts->report.unknown_barcodes++;
            continue;
        }
        else if (b->flag == FQ_PASS) {
            fprintf(fp1, "%c%s\n%s\n", b->q0 ? '@' : '>', b->n0, b->s0);
            if (b->q0) fprintf(fp1, "+\n%s\n", b->q0);
            if (b->l1 > 0) {
                fprintf(fp2, "%c%s\n%s\n", b->q1 ? '@' : '>', b->n0, b->s1);
                if (b->q1) fprintf(fp2, "+\n%s\n", b->q1);
            }
        
        }
        if (b->data) {
            if (opts->barcode_fp) {
                fputs((char*)b->data, opts->barcode_fp);
                fputc('\n', opts->barcode_fp);
            }
            free(b->data);
        }
    }
    
    bseq_pool_destroy(p);
}
int parse_args(int argc, char **argv)
{
    if (argc == 1) {
        fprintf(stderr, "* %s -platform [sci|drop] read_1.fq.gz read_2.fq.gz\n", program_name);
        fprintf(stderr, "* Use -h for more information.\n");
        return 1;
    }

    int i;
    const char *thread = NULL;
    const char *chunk_size = NULL;
    const char *platform = NULL;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return usage();
        if (strcmp(a, "-v") == 0 || strcmp(a, "--version") == 0) return get_version();
        if (strcmp(a, "-1") == 0) var = &args.out1_fname;
        else if (strcmp(a, "-2") == 0) var = &args.out2_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-platform") == 0) var = &platform;
        else if (strcmp(a, "-lane") == 0) var = &args.lane;
        else if (strcmp(a, "-cb") == 0) var = &args.barcode_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-r") == 0) var = &chunk_size;
        else if (strcmp(a, "-name") == 0) var = &args.name;
        else if (strcmp(a, "-q") == 0) {
            args.rm_qc_failed = 0;
            continue;
        }
        else if (strcmp(a, "-b") == 0) {
            args.rm_unknown_bc = 0;
            continue;
        }
        else if (strcmp(a, "-p") == 0) {
            args.smart_pair = 1;
            continue;
        }
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }
        if (a[0] == '-' && a[1]) error("Unknown parameter %s.", a);
        if (args.r1_fname == 0) {
            args.r1_fname = a;
            continue;
        }
        if (args.r2_fname == 0) {
            args.r2_fname = a;
            continue;
        }
        error("Unknown argument: %s, use -h see help information.", a);
    }

    if (thread) args.n_thread = str2int((char*)thread);
    if (chunk_size) args.chunk_size = str2int((char*)chunk_size);
    assert(args.n_thread >= 1 && args.chunk_size >= 1);

    if (args.name != NULL ) args.platform_id = PLT_RENAME;
    else {
        if (platform == NULL) error("-platform must be set.");
        else {
            if (strcmp(platform, "drop") == 0) args.platform_id = PLT_DROP;
            else if (strcmp(platform, "sci") == 0) args.platform_id = PLT_SCI;
            else if (strcmp(platform, "stream") == 0 ) args.platform_id = PLT_STREAM;
            else if (strcmp(platform, "plate") == 0 ) args.platform_id = PLT_PLATE;
            else error("Unknown platform: %s.", platform);
        }
    }
    if (args.r1_fname == NULL && (!isatty(fileno(stdin)))) args.r1_fname = "-";
    if (args.r1_fname == NULL) error("Fastq file(s) must be set.");

    if (args.out1_fname) {
        args.out1_fp = fopen(args.out1_fname, "w");
        if (args.out1_fp == NULL) error("%s: %s.", args.out1_fname, strerror(errno));
        if (args.out2_fname) {
            args.out2_fp = fopen(args.out2_fname,"w");
            if (args.out2_fp == NULL) error("%s: %s.", args.out2_fname, strerror(errno));
        }
    }

    args.fastq = fastq_handler_init(args.r1_fname, args.r2_fname, args.smart_pair, args.chunk_size);
    if (args.fastq == NULL) error("Failed to init input fastq.");
    
    if (args.report_fname) {
        args.report_fp = fopen(args.report_fname, "w");
        if (args.report_fp == NULL) error("%s: %s.", args.report_fname, strerror(errno));
    }
    if (args.barcode_fname) {
        args.barcode_fp = fopen(args.barcode_fname, "w");
        if (args.barcode_fp == NULL) error("%s: %s.", args.barcode_fname, strerror(errno));
    }
    return 0;
}
void write_report()
{
    if (args.report_fp) {
        fprintf(args.report_fp, "| Name | Value |\n|----|:----:|\n");
        fprintf(args.report_fp, "| All fragments | %llu |\n", args.report.all_fragments);
        fprintf(args.report_fp, "| QC failed fragments | %llu |\n", args.report.qc_failed);
        fprintf(args.report_fp, "| Unknown barcodes fragments | %llu |\n", args.report.unknown_barcodes);
        fclose(args.report_fp);
    }
}
void memory_release()
{
    if (args.r1_fp) gzclose(args.r1_fp);
    if (args.r2_fp) gzclose(args.r2_fp);
    if (args.out1_fp) fclose(args.out1_fp);
    if (args.out2_fp) fclose(args.out2_fp);
    //if (args.report_fp) fclose(args.report_fp);
    if (args.barcode_fp) fclose(args.barcode_fp);
    fastq_handler_destory(args.fastq);
}
int main(int argc, char **argv)
{
    if (parse_args(argc, argv)) return 1;

    void *opts = &args;
    
    if (args.platform_id == PLT_SCI)
        multi_threads_workflow(sci_run, write_out, (void*)args.fastq, bseq_read, opts, 0, args.n_thread);
    else if (args.platform_id == PLT_DROP)
        multi_threads_workflow(drop_run, write_out, (void*)args.fastq, bseq_read, opts, 0, args.n_thread);
    else if (args.platform_id == PLT_RENAME)
        multi_threads_workflow(manually_rename, write_out, (void*)args.fastq, bseq_read, opts, 0, args.n_thread);
    else if (args.platform_id == PLT_STREAM)
        multi_threads_workflow(stream_run, write_out, (void*)args.fastq, bseq_read, opts, 0, args.n_thread);
    else if (args.platform_id == PLT_PLATE)
        multi_threads_workflow(plate_run, write_out, (void*)args.fastq, bseq_read, opts, 0, args.n_thread);
    write_report();
    
    memory_release();
    return 0;
}
