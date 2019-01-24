#include "utils.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"

int usage()
{
    fprintf(stderr, "cb_count_matrix [options] aln.bam\n");
    fprintf(stderr, "options:\n");
    fprintf(stderr, "  -report report.txt     Summary report.\n");
    fprintf(stderr, "  -barcode barcode.txt   Predefined barcode list, usually be defined by QC steps.\n");
    fprintf(stderr, "  -tag CR                Cell barcode tag at BAM file.\n");
    fprintf(stderr, "  -bed target.bed        Peak file or gene region.\n");
    fprintf(stderr, "  -n                     Use peak name as row name, merge sub-regions with same name.\n");
    fprintf(stderr, "  -t 5                   Thread.\n");
    return 1;
}
struct args {
    const char *input_fname;
    const char *report_fname;
    const char *stat_fname;
} args = {
};
struct barcode_str {
    char *str;
    int count;
};

void cell_barcode_stat()
{
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) ) return usage();

    return cell_barcode_stat();
}
