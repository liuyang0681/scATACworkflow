task align {
  String Root
  String reference
  String Outdir
  String FQ1
  String FQ2
  String name
  command {
    ${Root}/third_party/bwa mem -t 20  ${reference} ${FQ1} ${FQ2} | \
    ${Root}/bin/parse_alignment_readname -filter ${Outdir}/temp/${name}_filter.bam -o ${Outdir}/temp/${name}_aln.bam -report ${Outdir}/temp/${name}_aln.md -t 5 -maln ${Outdir}/temp/${name}_mito.bam
  }
  output {
    String bamfile="${Outdir}/temp/${name}_aln.bam"
    String aln_report="${Outdir}/temp/${name}_aln.md"
    String mito_aln="${Outdir}/temp/${name}_mito.bam"
  }
}

workflow alignment {
  String Root
  String reference
  String Outdir
  String FQ1
  String FQ2
  String name
  call align {
    input:
    Root=Root,
    reference=reference,
    Outdir=Outdir,
    FQ1=FQ1,
    FQ2=FQ2,
    name=name,
  }
}
