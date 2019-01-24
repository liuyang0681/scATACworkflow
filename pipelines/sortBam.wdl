task sortBAM {
  String Outdir
  String Root
  String name
  String bamfile  
  command {
    ${Root}/third_party/sambamba sort -t 24 -o ${Outdir}/temp/${name}_sorted.bam ${bamfile}
    java -jar ${Root}/third_party/picard.jar MarkDuplicates I=${Outdir}/temp/${name}_sorted.bam O=${Outdir}/alignment/${name}_rmdup.bam M=${Outdir}/temp/${name}_dup.matrix Barcode_tag=CR REMOVE_DUPLICATES=true
  }
  output {
    String final="${Outdir}/alignment/${name}_rmdup.bam"
  }
}

