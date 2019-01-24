//
//
//
//
//
#include "utils.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include <zlib.h>
#include "thread_pool.h"
#include "number.h"

KSTREAM_INIT(gzFile, gzread, 16384)

// flag to skip
#define FLG_USABLE 0
#define FLG_MITO 1
#define FLG_FLT  2

static char *version = "0.0.0.9000";

static int get_version()
{
    puts(version);
    return 1;
}

// summary structure, thread safe

// cell barcode counts
struct cr_count {    
    char *str;
    int count;    
};
// count pool
struct cr {
    int n, m;
    struct cr_count *cr;
};

void insert_cr(struct cr *r, int start, int end, char *str, int c)
{
    int i;
    for (i = start; i <= end; ++i) {
        int v = strcmp(r->cr[i].str, str);
        if (v == 0) {
            r->cr[i].count += c;
            return;
        }
        else if (v > 0) {
            memmove(r->cr+i+1, r->cr+i, (r->n-i)*sizeof(struct cr_count));
            r->cr[i].str = strdup(str);
            r->cr[i].count = c;
            r->n++;
            return;
        }
    }    
    r->cr[r->n].str = strdup(str);
    r->cr[r->n].count = c;
    r->n++;
    
    /*
    assert(start <= end);
    
    int i;
    i =  strcmp(r->cr[start].str, str);
    if (i == 0) { // found at start
        r->cr[start].count += c;
        return;
    }
    else if (i < 0) goto extern_before_start;

    if (start == end) goto extern_at_end;

    //if (end == start +1) goto extern_before_start;

    i = strcmp(r->cr[end].str, str);
    if (i == 0) { // found at end
        r->cr[end].count += c;
        return;
    }
    else if (i > 0) goto extern_at_end;
    
    int m = (start+end)/2+1;
    i = strcmp(r->cr[m].str, str);
    if (i == 0) {
        r->cr[m].count += c;
        return;
    }
    else if (i < 0) {
        return insert_cr(r, start+1, m, str, c);
    }
    else {
        return insert_cr(r, m+1, end, str, c);
    }

    
  extern_at_end:
    memmove(r->cr+end+1, r->cr+end, (r->n - end)*sizeof(struct cr_count));
    r->cr[end].str = strdup(str);
    r->cr[end].count = c;
    r->n++;
    return;

  extern_before_start:
    memmove(r->cr+start+1, r->cr+start, (r->n - start)*sizeof(struct cr_count));
    r->cr[start].str = strdup(str);
    r->cr[start].count = c;
    r->n++;
    return;
    */
}
void push_cr(struct cr *r, char *str, int c)
{
    // debug_print("%p\t%s\t%d\t%d\t%d", r, str, c, r->n, r->m);
    if (r->n == r->m) {
        r->m = r->m == 0 ? 1024 : r->m*2;
        r->cr = realloc(r->cr, r->m *sizeof(struct cr_count));
    }
        
    if (r->n == 0) {
        r->cr[0].str = strdup(str);
        r->cr[0].count = c;
        r->n++;
    }
    else
        insert_cr(r, 0, r->n-1, str, c);
}
// add barcode information cr1 to cr
void add_cr(struct cr *cr, struct cr *cr1)
{
    int i;
    for (i = 0; i < cr1->n; ++i) {
        struct cr_count *c = &cr1->cr[i];
        push_cr(cr, c->str, c->count);
    }
}

// summary structure for final report
struct reads_summary {
    uint64_t n_reads, n_mapped, n_pair_map, n_pair_all, n_pair_good;
    uint64_t n_sgltn, n_read1, n_read2;
    uint64_t n_diffchr, n_pstrand, n_mstrand;
    uint64_t n_qual;
    uint64_t n_mito;
    uint64_t n_usable;
    struct cr *cr;
};
struct reads_summary *reads_summary_create()
{
    struct reads_summary *s = malloc(sizeof(*s));
    memset(s, 0, sizeof(*s));
    s->cr = malloc(sizeof(struct cr));
    memset(s->cr, 0, sizeof(struct cr));
    return s;
}
struct reads_summary *merge_reads_summary(struct reads_summary **rs, int n)
{
    struct reads_summary *s = reads_summary_create();

    int i;
    for (i = 0; i < n; ++i) {

        struct reads_summary *s0 = rs[i];

        s->n_reads     += s0->n_reads;
        s->n_mapped    += s0->n_mapped;
        s->n_pair_map  += s0->n_pair_map;
        s->n_pair_all  += s0->n_pair_all;
        s->n_pair_good += s0->n_pair_good;
        s->n_sgltn     += s0->n_sgltn;
        s->n_read1     += s0->n_read1;
        s->n_read2     += s0->n_read2;
        s->n_diffchr   += s0->n_diffchr;
        s->n_pstrand   += s0->n_pstrand;
        s->n_mstrand   += s0->n_mstrand;
        s->n_qual      += s0->n_qual;
        s->n_mito      += s0->n_mito;
        s->n_usable    += s0->n_usable;
        add_cr(s->cr, s0->cr);
    }
    return s;
}
void reads_summary_destory(struct reads_summary *s)
{
    int i;
    struct cr *cr = s->cr;
    for (i = 0; i < cr->n; ++i) 
        free(cr->cr[i].str);
    free(cr->cr);
    free(cr);
    free(s);
}

struct args {
    // file names
    const char *input_fname;  // alignment input, only SAM format is required
    const char *output_fname; // BAM output, only support BAM format, default is stdout
    const char *filter_fname; // if filter BAM file is set, filtered reads will output in this file
    const char *report_fname; // summary report for whole file
    //const char *btable_fname; // export raw read count for each barcode, sorting according to read counts

    const char *mito; // mitochrondria name, the mito ratio will export in the summary report
    const char *mito_fname; // if set, mitochrondria reads passed QC will be exported in this file
    
    int n_thread;
    int buffer_size;  // buffered records in each chunk
    
    gzFile fp;        // input file handler
    kstream_t *ks;    // input streaming
    BGZF *fp_out;     // output file handler
    BGZF *fp_filter;  // filtered reads file handler
    BGZF *fp_mito;    // if not set, mito reads will be treat at filtered reads
    FILE *fp_report;  // report file handler
    //FILE *fp_btable;  // barcode list table 
    
    bam_hdr_t *hdr;   // bam header structure of input

    struct reads_summary **thread_data; // summary data for each thread
    
    // struct reads_summary *summary; // final summary for report and table

    int fixmate_flag; // if set, fix the mate relationship of paired reads first

    int qual_thres; // mapping quality threshold to filter

    int mito_id;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .filter_fname = NULL,
    .report_fname = NULL,
    //.btable_fname = NULL,
    .mito = "chrM",
    .mito_fname = NULL,
    
    .n_thread = 5,
    .buffer_size = 1000000, // 1M

    .fp = NULL,
    .ks = NULL,
    .fp_out = NULL,
    .fp_filter = NULL,
    .fp_mito = NULL,
    .fp_report = NULL,
    //.fp_btable = NULL,
    
    .hdr = NULL,

    .thread_data = NULL,
    //.summary = NULL,

    .fixmate_flag = 0,
    .qual_thres = 10,
    .mito_id = -2,
};

// Buffer input and output records in a memory pool per thread
struct sam_pool {
    struct args *opts; // point to args
    int n;
    kstring_t **str;
    bam1_t **bam; // bam structure
    int *flag; // export flag
};

struct sam_pool* sam_pool_init(int buffer_size)
{
    struct sam_pool *p = malloc(sizeof(*p));
    p->n = 0;

    p->str  = malloc(buffer_size*sizeof(void*));
    p->bam  = malloc(buffer_size*sizeof(void*));
    p->flag = malloc(buffer_size*sizeof(int));
    
    memset(p->str,  0, buffer_size*sizeof(void*));
    memset(p->bam,  0, buffer_size*sizeof(void*));
    memset(p->flag, 0, buffer_size*sizeof(int));
    
    return p;
}
void sam_pool_destroy(struct sam_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i) {
        if (p->str[i]) {
            if (p->str[i]->m) free(p->str[i]->s);
            free(p->str[i]);
        }
        if (p->bam[i]) bam_destroy1(p->bam[i]);        
    }
    free(p->str);
    free(p->bam);
    free(p->flag);
    free(p);
}
struct sam_pool* sam_pool_read(kstream_t *s, int buffer_size)
{
    struct sam_pool *p = sam_pool_init(buffer_size);
    
    kstring_t str = {0,0,0};
    int ret;
    
    for (;;) {
        
        if (ks_getuntil(s, 2, &str, &ret) < 0) break;
        
        // skip header
        if (str.s[0] == '@') continue;

        kstring_t *t = malloc(sizeof(*t));
        t->l = str.l;
        t->m = str.m;
        t->s = strdup(str.s);
        p->str[p->n] = t;
        // init for output
        p->bam[p->n] = bam_init1();
        
        p->n++;
        if (p->n >= buffer_size) break;
    }

    free(str.s);
    if (p->n == 0) {
        sam_pool_destroy(p);
        return NULL;
    }
    
    return p;
}
bam_hdr_t *sam_parse_header(kstream_t *s, kstring_t *line)
{
    bam_hdr_t *h = NULL;
    kstring_t str = {0,0,0};
    int ret;
    
    while (ks_getuntil(s, 2, line, &ret) >= 0) {
        if (line->s[0] != '@') break;
        kputsn(line->s, line->l, &str);
        kputc('\n', &str);
    }

    if (ret < -1) goto failed_parse;
    if (str.l == 0) kputsn("", 0, &str);
    h = sam_hdr_parse(str.l, str.s);
    h->l_text = str.l;
    h->text = str.s;
    return h;

  failed_parse:
    bam_hdr_destroy(h);
    free(str.s);
    
    return NULL;
}
void write_out(struct sam_pool *p)
{
    int i;
    struct args *opts = p->opts;
    for (i = 0; i < p->n; ++i) {
        if (p->flag[i] == FLG_USABLE ) {
            if (bam_write1(opts->fp_out, p->bam[i]) == -1) error("Failed to write.");            
        }
        else if (p->flag[i] == FLG_MITO && opts->fp_mito != NULL) {
            if (bam_write1(opts->fp_mito, p->bam[i]) == -1) error("Failed to write.");
        }
        else if (opts->fp_filter != NULL) {
            if (bam_write1(opts->fp_filter, p->bam[i]) == -1) error("Failed to write.");
        }
    }
    sam_pool_destroy(p);
}
void summary_report(struct args *opts)
{
    struct reads_summary *summary = merge_reads_summary(opts->thread_data, opts->n_thread);
    if (opts->fp_report) {
        fprintf(opts->fp_report,"| Name | Value |\n|----|:----:|\n");
        fprintf(opts->fp_report,"| Raw reads | %"PRIu64" |\n", summary->n_reads);
        fprintf(opts->fp_report,"| Mapped reads | %"PRIu64" (%.2f%%) |\n", summary->n_mapped, (float)summary->n_mapped/summary->n_reads*100);
        fprintf(opts->fp_report,"| Mapped reads (paired) | %"PRIu64" |\n", summary->n_pair_map);
        fprintf(opts->fp_report,"| Properly paired reads | %"PRIu64" |\n", summary->n_pair_good);
        fprintf(opts->fp_report,"| Singleton reads | %"PRIu64" |\n", summary->n_sgltn);
        fprintf(opts->fp_report,"| Read 1 | %"PRIu64" |\n", summary->n_read1);
        fprintf(opts->fp_report,"| Read 2 | %"PRIu64" |\n", summary->n_read2);
        fprintf(opts->fp_report,"| Paired reads map on diff chr | %"PRIu64" |\n", summary->n_diffchr);
        fprintf(opts->fp_report,"| Plus strand | %"PRIu64" |\n", summary->n_pstrand);
        fprintf(opts->fp_report,"| Minus strand | %"PRIu64" |\n", summary->n_mstrand);
        fprintf(opts->fp_report,"| Mapping quals above %d | %"PRIu64" |\n", opts->qual_thres, summary->n_qual);
        fprintf(opts->fp_report,"| Mitochondria ratio | %.2f%% |\n", (float)summary->n_mito/summary->n_mapped*100);
        fprintf(opts->fp_report,"| Usable reads (ratio) | %"PRIu64" (%.2f%%)|\n", summary->n_usable, (float)summary->n_usable/summary->n_reads*100);
    }
    /*
    if (opts->fp_btable) {
        int i;
        struct cr *cr = summary->cr;
        for (i = 0; i < cr->n; ++i) {
            fprintf(opts->fp_btable, "%s\t%d\n", cr->cr[i].str, cr->cr[i].count);
        }
    }
    */
    reads_summary_destory(summary);
}

static char* parse_name_str(kstring_t *s)
{
    // CL100053545L1C001R001_2|||BC|||TTTCATGA|||CR|||TANTGGTAGCCACTAT
    // CL100053545L1C001R001_2 .. CR:Z:TANTGGTAGCCACTAT ..
    int n, i;
    for (n = 0; n < s->l && !isspace(s->s[n]); ++n);
    for (i = 0; i < n && s->s[i] != '|'; ++i);
    
    int is_name = 1; // first key is name, next is value
    char *cr_name = NULL;
    if (i < n-8) { // detected
        char *p = s->s+i;
        char *e = s->s+n;
        kstring_t t = {0,0,0};
        kputsn(s->s, i, &t);
        kputsn(s->s+n, s->l-n, &t);
        int is_cr = 0;
        for (;;) {
            if (*p == '|' && *(p+1) == '|' && *(p+2) == '|') {
                p = p + 3;
                if (is_name) {
                    if (*(p + 2) != '|') {
                        warnings("Failed to parse name %s", s->s);
                        if (t.m) free(t.s);
                        break;
                    }
                    kputc('\t', &t); kputsn(p, 2, &t); kputs(":Z:", &t);

                    if (strncmp(p, "CR", 2) == 0) is_cr = 1;
                    
                    p = p + 2;
                    is_name = 0;
                }
                else {
                    char *pp = p;
                    int l = 0;
                    while (*pp != '|' && pp != e) { ++pp; ++l; }
                    kputsn(p, l, &t);

                    if (is_cr) {
                        cr_name = strndup(p, l);
                        is_cr = 0;
                    }                        

                    is_name = 1;
                    p = pp;
                }
                continue;
            }
            break;
        }
        s->l = t.l;
        s->m = t.m;
        free(s->s);
        s->s = t.s;
    }
    return cr_name;
}
void count_cell_barcode(char *name, struct reads_summary *s)
{
    struct cr *cr = s->cr;
    push_cr(cr, name, 1); // raw reads, do not apply any filter
}
void sam_stat_reads(bam1_t *b, struct reads_summary *s, int *flag, struct args *opts)
{
    bam1_core_t *c = &b->core;
    s->n_reads++;
    if (c->flag & BAM_FQCFAIL || c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY || c->flag & BAM_FUNMAP)
        *flag = FLG_FLT; // filtered
    else {
        // mito ratio, mito come before qc filter, so only high quality mito reads will be exported in mito.bam
        if (c->tid == opts->mito_id) {
            s->n_mito++;
            *flag = FLG_MITO; // filter mito
        }
        
        s->n_mapped++;
        if (c->flag & BAM_FREVERSE) s->n_mstrand++;
        else s->n_pstrand++;            

        if (c->qual >= opts->qual_thres) s->n_qual++;
        else *flag = FLG_FLT; // filter low mapping quality
        
        if (c->flag & BAM_FPAIRED) {
            s->n_pair_all ++;
            if ((c->flag & BAM_FPROPER_PAIR) && !(c->flag & BAM_FUNMAP)) s->n_pair_good++;
            if (c->flag & BAM_FREAD1) s->n_read1++;
            else if (c->flag & BAM_FREAD2) s->n_read2++;
            if ((c->flag & BAM_FMUNMAP) && !(c->flag & BAM_FUNMAP)) s->n_sgltn++;
            if (!(c->flag & BAM_FUNMAP) && !(c->flag & BAM_FMUNMAP)) {
                s->n_pair_map++;
                if (c->mtid != c->tid) s->n_diffchr++;
            }
        }
    }
    if (*flag == FLG_USABLE) s->n_usable++;
}
void *sam_name_parse(void *_p, int idx)
{
    struct sam_pool *p = (struct sam_pool*)_p;
    struct args *opts = p->opts;
    struct reads_summary *summary = opts->thread_data[idx];    
    bam_hdr_t *h = opts->hdr;
    
    int i;
    for (i = 0; i < p->n; ++i) {
        parse_name_str(p->str[i]);
        //count_cell_barcode(name, summary);
        //free(name);
        if (sam_parse1(p->str[i], h, p->bam[i])) error("Failed to parse SAM.");
        sam_stat_reads(p->bam[i], summary, &p->flag[i], opts);
    }
    
    return p;
}
int sam_name_parse_light()
{
    for (;;) {
        struct sam_pool *p = sam_pool_read(args.ks, args.buffer_size);
        if (p == NULL) break;
        p->opts = &args;
        p = sam_name_parse(p, 0);
        write_out(p);
    }
    return 0;
}
int usage()
{
    fprintf(stderr, "name_parser in.sam\n");
    fprintf(stderr, " -o out.bam               Output file [stdout].\n");
    fprintf(stderr, " -t [5]                   Threads.\n");
    fprintf(stderr, " -filter filter.bam       Filter reads, include unmapped, secondary alignment, mitochondra etc.\n");
    fprintf(stderr, " -q [10]                  Map quality theshold, mapped reads with smaller mapQ will be filter.\n");
    fprintf(stderr, " -report report.md        Summary report in Markdown format. Alignment stats, chromosome ratio.\n");
    //fprintf(stderr, " -bsum bar_sum.txt        Summary report for barcode. Export raw reads count for each barcode.\n");
    fprintf(stderr, " -mito chrM               Mitochondria name.\n");
    fprintf(stderr, " -maln chrM.bam           Export mitochondria reads into this file.\n");
    fprintf(stderr, " -r [1000000]             Records per chunk.\n");
    return 1;
}

int bam_fixmate()
{
    return 0;
}

int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;
    
    int i;
    const char *buffer_size = NULL;
    const char *thread = NULL;
    const char *qual_thres = NULL;
    
    for (i = 1; i < argc;) {

        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0)
            return usage();
        if (strcmp(a, "-v") == 0 || strcmp(a, "--version") == 0 )
            return get_version();
        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-r") == 0) var = &buffer_size;
        else if (strcmp(a, "-q") == 0) var = &qual_thres;
        else if (strcmp(a, "-mito") == 0) var = &args.mito;
        else if (strcmp(a, "-maln") == 0) var = &args.mito_fname;
        // else if (strcmp(a, "-bsum") == 0) var = &args.btable_fname;
        else if (strcmp(a, "-report") == 0) var = &args.report_fname;
        else if (strcmp(a, "-filter") == 0) var = &args.filter_fname;

        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1]) error("Unknown parameter: %s", a);

        if (args.input_fname == 0) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument: %s.", a);
    }

    // init input
    if (args.input_fname == NULL && !isatty(fileno(stdin))) args.input_fname = "-";
    if (args.input_fname == NULL) error("No input SAM file is set!");
    args.fp = strcmp(args.input_fname, "-") ? gzopen(args.input_fname, "r") : gzdopen(fileno(stdin), "r");
    if (args.fp == NULL) error("%s : %s.", args.input_fname, strerror(errno));
    args.ks = ks_init(args.fp);

    // init output
    args.fp_out = args.output_fname ? bgzf_open(args.output_fname, "w") : bgzf_fdopen(fileno(stdout), "w");
    if (args.fp_out == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    /*
    if (args.btable_fname) {
        args.fp_btable = fopen(args.btable_fname, "w");
        if (args.fp_btable == NULL) error("%s : %s.", args.btable_fname, strerror(errno));
    }
    */
    if (args.report_fname) {
        args.fp_report = fopen(args.report_fname, "w");
        if (args.fp_report == NULL) error("%s : %s.", args.report_fname, strerror(errno));
    }
    if (args.filter_fname) {
        args.fp_filter = bgzf_open(args.filter_fname, "w");
        if (args.fp_filter == NULL) error("%s : %s.", args.filter_fname, strerror(errno));
    }
    if (args.mito_fname) {
        args.fp_mito = bgzf_open(args.mito_fname, "w");
        if (args.fp_mito == NULL) error("%s : %s.", args.mito_fname, strerror(errno));
    }
    // init parameters
    if (thread) {
        args.n_thread = str2int((char*)thread);
        assert(args.n_thread > 0);
    }
    if (buffer_size) {
        args.buffer_size = str2int((char*)buffer_size);
        assert(args.buffer_size>0);
    }
    if (qual_thres) {
        args.qual_thres = str2int((char*)qual_thres);
        assert(args.qual_thres >= 0);
    }

    // init thread data
    args.thread_data = malloc(args.n_thread*sizeof(struct reads_summary*));
    for (i = 0; i < args.n_thread; ++i) args.thread_data[i] = reads_summary_create();

        
    // init bam header and first bam record
    kstring_t str = {0,0,0}; // cache first record
    args.hdr = sam_parse_header(args.ks, &str);
    if (args.hdr == NULL) error("Failed to parse header. %s", args.input_fname);
    if (bam_hdr_write(args.fp_out, args.hdr)) error("Failed to write header.");
    if (args.fp_filter && bam_hdr_write(args.fp_filter, args.hdr)) error("Failed to write header.");
    if (args.fp_mito && bam_hdr_write(args.fp_mito, args.hdr)) error("Failed to write header.");

    // init mitochondria id
    args.mito_id = bam_name2id(args.hdr, args.mito);
    if (args.mito_id == -1) {
        warnings("Failed to find mitochondria name at the header, will NOT stat the mito ratio!");
        args.mito_id = -2; // reset to -2, in case inconsistance with unmap
    }

    // check if there is a BAM record
    if (str.s[0] != '@') {
        parse_name_str(&str);
        bam1_t *b = bam_init1();
        if (sam_parse1(&str, args.hdr, b)) error("Failed to parse SAM.");

        /*
        if (cr_name) {
            // add this record to first thread data
            struct cr *cr = args.thread_data[0]->cr;
            push_cr(cr, cr_name, 1);
            free(cr_name);
        }
        */
        int flag = 0;        
        sam_stat_reads(b, args.thread_data[0], &flag, &args);
        
        if (flag == FLG_USABLE) {
            if (bam_write1(args.fp_out, b) == -1) error("Failed to write BAM.");
        }
        else if (flag == FLG_MITO && args.fp_mito != NULL) {
            if (bam_write1(args.fp_mito, b) == -1) error("Failed to write BAM.");
        }
        else if (args.fp_filter != NULL) {
            if (bam_write1(args.fp_filter, b) == -1) error("Failed to write BAM.");
        }
        else {

        }
        bam_destroy1(b);
    }
    free(str.s);
    return 0;
}

void memory_release()
{
    bgzf_close(args.fp_out);
    ks_destroy(args.ks);
    gzclose(args.fp);
    bam_hdr_destroy(args.hdr);

    if (args.fp_filter) bgzf_close(args.fp_filter);
    if (args.fp_mito) bgzf_close(args.fp_mito);
    if (args.fp_report) fclose(args.fp_report);
    //if (args.fp_btable) fclose(args.fp_btable);

    int i;
    for (i = 0; i < args.n_thread; ++i) reads_summary_destory(args.thread_data[i]);
    free(args.thread_data);
}

int main(int argc, char **argv)
{
    double t_real;
    t_real = realtime();
    
    if (parse_args(argc, argv)) return usage();

    if (args.n_thread == 1)
        sam_name_parse_light();
    
    else {

        int nt = args.n_thread;

        struct thread_pool *p = thread_pool_init(nt);
        struct thread_pool_process *q = thread_pool_process_init(p, nt*2, 0);
        struct thread_pool_result *r;

        for (;;) {

            struct sam_pool *b = sam_pool_read(args.ks, args.buffer_size);
            if (b == NULL) break;
            b->opts = &args;
            
            int block;

            do {

                block = thread_pool_dispatch2(p, q, sam_name_parse, b, 0);

                if ((r = thread_pool_next_result(q))) {
                    struct sam_pool *d = (struct sam_pool*)r->data;
                    write_out(d);
                }
                //thread_pool_delete_result(r, 1);
            }
            while (block == -1);
        }
        
        thread_pool_process_flush(q);

        while ((r = thread_pool_next_result(q))) {
            struct sam_pool *d = (struct sam_pool *)r->data;
            write_out(d);                
        }

        thread_pool_process_destroy(q);
        thread_pool_destroy(p);
        
    }

    summary_report(&args);
    
    memory_release();

    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}
