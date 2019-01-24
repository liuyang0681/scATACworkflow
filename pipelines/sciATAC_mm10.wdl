import "scATAC_core.wdl" as core

workflow main {
  String Root
  String Fastq1
  String Fastq2
  String Outdir
  String TestID
  String? lane
  call core.core {
    input:
    Root=Root,
    FQ1=Fastq1,
    FQ2=Fastq2,
    Outdir=Outdir,
    mode="sci",
    reference="${Root}/reference/mm10/mm10_no_alt_analysis_set_ENCODE.fasta",
    tssflank="${Root}/file/mm10_gencode_tss_flank_r2K+l2K.bed",
    tssbed="${Root}/file/mm10_gencode_tss_unique.bed",
    gsize="mm",
    TestID=TestID,
    lane=lane,
  }
}
