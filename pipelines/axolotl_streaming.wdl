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
    mode="stream",
    reference="${Root}/reference/axolotl/axolotl.fa",
    tssflank="${Root}/file/axolotl_tss_flank_r2K+l2K.bed",
    tssbed="${Root}/file/axolotl_tss_unique.bed",
    gsize="32406177977",
    TestID=TestID,
    lane=lane,
  }
}
