task accessEnhancer {
  String Root
  String TssBed
  String name
  String finalBAM
  String Outdir
  command <<<
    ${Root}/third_party/bamCoverage -b ${finalBAM} -o ${Outdir}/report/${name}.bw -of bigwig &&\
    ${Root}/third_party/computeMatrix reference-point -S ${Outdir}/report/${name}.bw \
    -R ${TssBed} -b 3000 -a 3000 --skipZeros -o ${Outdir}/temp/${name}_tss_access.matrix && \
    ${Root}/third_party/plotHeatmap -m ${Outdir}/temp/${name}_tss_access.matrix \
    -out ${Outdir}/temp/${name}_tss_access.pdf --colorMap GnBu --plotTitle '${name} cells'
  >>>
  output {
    String accEnhPlot="${Outdir}/temp/${name}_tss_access.pdf"
  }
}
