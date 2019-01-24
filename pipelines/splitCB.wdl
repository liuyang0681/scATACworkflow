task split {
  String Outdir
  String Fastq1
  String Fastq2
  String Root
  String mode
  command {
    ${Root}/bin/sc_atac_parse_cellbarcode -t 10 ${Fastq1} ${Fastq2} -platform ${mode} -report ${Outdir}/temp/cb_report.md -cb ${Outdir}/temp/barcode_table.txt| \
    ${Root}/bin/parse_Tn5me_adaptors -t 10 -d -p -report ${Outdir}/temp/trim_me_report.md -1 ${Outdir}/fastq/read_1.fq -2 ${Outdir}/fastq/read_2.fq
  }
  output {
    String CBtable="${Outdir}/temp/barcode_table.txt"
    String FQ1="${Outdir}/fastq/read_1.fq"
    String FQ2="${Outdir}/fastq/read_2.fq"
    String TMreport="${Outdir}/temp/trim_me_report.md"
    String CBreport="${Outdir}/temp/cb_report.md"
  }
}

workflow splitCB {
  String Outdir
  String Fastq1
  String Fastq2
  String Root
  call split {
    input:
    Outdir=Outdir,
    Fastq1=Fastq1,
    Fastq2=Fastq2,
    Root=Root,
  }

  output {
    String CBtable=split.CBtable
    String FQ1=split.FQ1
    String FQ2=split.FQ2
  }
}
