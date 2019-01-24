task report {
  String Outdir
  String Root
  String matrix
  String Sample
  String accEnhPlot
  String frag
  String CBreport
  String TMreport
  String ALNreport
  command <<<
        echo '---
title: "scATAC pipeline, ${Sample}"
output:
  html_document:
    toc: true
    theme: united
  author: "scATAC-pipeline automatical reportor"
  date: "`r format(Sys.time(), '"'"'%d %B, %Y'"'"')`"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(root.dir = normalizePath("${Outdir}/report"))
```
' > ${Outdir}/report/report.Rmd

    echo "" >>${Outdir}/report/report.Rmd
    echo "# Split cell barcode" >> ${Outdir}/report/report.Rmd
    cat ${CBreport} >> ${Outdir}/report/report.Rmd

    echo "" >>${Outdir}/report/report.Rmd
    echo "# Trim Tn5 Mosaic Ends and adaptors" >> ${Outdir}/report/report.Rmd
    cat ${TMreport} >> ${Outdir}/report/report.Rmd

    echo "" >>${Outdir}/report/report.Rmd
    echo "# Alignment state" >>${Outdir}/report/report.Rmd   
    cat ${ALNreport} >> ${Outdir}/report/report.Rmd

    echo "" >>${Outdir}/report/report.Rmd

echo "# Fragement size distribution" >> >>${Outdir}/report/report.Rmd
    echo '```{r echo=FALSE, fig.height = 3,  fig.width = 6, fig.align = "center"}
library("ggplot2")
library("gridExtra")
library("data.table")
  

    ```{r echo=FALSE, fig.height = 3,  fig.width = 12, fig.align = "center"}
    fg <- fread("${frag}", header=F)

    ggplot(data=fg, aes(V1, V2)) + geom_line()

    ```    
    ' >> ${Outdir}/report/report.Rmd
    
    export RSTUDIO_PANDOC=${Root}/third_party/pandoc
    ${Root}/third_party/R -e 'rmarkdown::render("${Outdir}/report/report.Rmd", output_format="html_document",output_file ="${Outdir}/report/report.html")'
    
    ${Root}/third_party/Rscript ${Root}/scripts/run_chromVAR.R ${matrix} ${Outdir}/report/chromVAR/${Sample}
  >>>
}
