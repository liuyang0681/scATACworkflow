task callpeak {
  String Outdir
  String Root
  String Sample
  String finalBAM
  String gsize
  command <<<
    ${Root}/third_party/macs2 callpeak -t ${finalBAM} -f BAM -n ${Sample} -q 0.01 --nomodel --shift -75 --extsize 150 -g ${gsize} --outdir ${Outdir}/temp/

    cut -f1,2,3 ${Outdir}/temp/${Sample}_peaks.narrowPeak | grep -v "_" > ${Outdir}/report/${Sample}_peak.bed
  >>>
  output {
    String peak="${Outdir}/report/${Sample}_peak.bed"
  }
}
