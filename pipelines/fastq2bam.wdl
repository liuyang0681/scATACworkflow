task fastq2bam {
  String Outdir
  String Fastq1
  String Fastq2
  String mode
  String Sample
  String reference
  String Root
  String? lane
  command {
    ${Root}/bin/sc_atac_parse_cellbarcode -t 10 ${Fastq1} ${Fastq2} -platform ${mode} -report ${Outdir}/temp/cb_report.md -cb ${Outdir}/temp/barcode_table.txt -lane ${default="1" lane}| \
    ${Root}/bin/parse_Tn5me_adaptors -t 10 -d -p -report ${Outdir}/temp/trim_me_report.md |\
    ${Root}/third_party/bwa mem -t 20 -p ${reference} /dev/stdin | \
    ${Root}/bin/parse_alignment_readname -filter ${Outdir}/temp/filter.bam -o ${Outdir}/temp/aln.bam -report ${Outdir}/temp/aln.md -t 5 -maln ${Outdir}/temp/mito.bam
  }
  output {
    String CBtable="${Outdir}/temp/barcode_table.txt"
    String mitoBAM="${Outdir}/temp/mito.bam"
    String rawBAM="${Outdir}/temp/aln.bam"
    String CBreport="${Outdir}/temp/cb_report.md"
    String TMreport="${Outdir}/temp/trim_me_report.md"
    String ALNreport="${Outdir}/temp/aln.md"
  }

}
