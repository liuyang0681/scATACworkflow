task countMTX {
  String Outdir
  String BClist
  String Root
  String peak
  command <<<
    python ${Root}/scripts/sc_atac_count.py ${Outdir}/alignment/aln_rmdup.bam ${BClist} ${peak} ${Outdir}/report/raw_peak_count_matrix.txt True
  >>>
  output {
    String matrix="${Outdir}/report/raw_peak_count_matrix.txt"
  }
}
