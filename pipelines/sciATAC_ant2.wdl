import "scATAC_core.wdl" as core

workflow main {
  String Root
  String Fastq1
  String Fastq2
  String Outdir
  String TestID  
  call core.core {
    input:
    Root=Root,
    FQ1=Fastq1,
    FQ2=Fastq2,
    Outdir=Outdir,
    mode="stream",
    reference="${Root}/reference/ant/ant.fa",
    tssflank="${Root}/file/ant_flank.bed",
    tssbed="${Root}/file/ant.bed",
    gsize="3.253727e+08",
    TestID=TestID,
  }
}
