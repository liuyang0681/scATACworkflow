import "hs_mm_mixed_core.wdl" as core

workflow main {
  String Root
  String TestID
  String Fastq1
  String Fastq2
  String Outdir
  call core.hs_mm_mixed_core {
    input:
    Root=Root,
    TestID=TestID,
    Fastq1=Fastq1,
    Fastq2=Fastq2,
    Outdir=Outdir,
    mode="drop"
  }  
}
