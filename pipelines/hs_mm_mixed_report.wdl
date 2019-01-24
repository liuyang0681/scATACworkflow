task report {
  String Root
  String Outdir
  String RawCount
  String CBreport
  String TMreport
  String HsAln
  String MmAln
  String mm_UsableCount
  String hs_UsableCount 
  String TestID
  String accEnh_hs_plot
  String accEnh_mm_plot
  String fragHs
  String fragMm
  command <<<
    perl -we '
    open RAW,"${RawCount}" or die;
    open MM,"${mm_UsableCount}" or die;
    open HS,"${hs_UsableCount}" or die;
    my %hash=();
    while(<RAW>){
      chomp;
      my ($id, $count)=split;
      $count = $count*2; # convert fragements to reads
      my @val=($count, 0,0,0,0);
      $hash{$id}=\@val;
    }
    close RAW;
    while(<HS>){
      chomp;
      my ($id, $us, $tss) = split;
      if (exists $hash{$id}) {
        $hash{$id}[1] = $us; $hash{$id}[2] = $tss;
      }
      else {
        print STDERR "$id at hs_rmdup.bam not shown in raw data.\n";
      }
    }
    close HS;
    while(<MM>){
      chomp;
      my ($id, $us, $tss) = split;
      if (exists $hash{$id}) {
        $hash{$id}[3]=$us; $hash{$id}[4] = $tss;
      }
      else {
        print STDERR "$id at mm_rmdup.bam not shown in raw data.\n";
      }
    }
    close MM;
    open OUT,">${Outdir}/report/result_${TestID}.txt" or die;
    print OUT "barcode_ID\tbarcode_reads\tHs_usable_reads\tHs_reads_in_promoters\tMm_usable_reads\tMm_reads_in_promoters\n";
    for my $key (keys %hash) {
      if ($hash{$key}[1] > 100 || $hash{$key}[3]>100) {
        print OUT "$key\t$hash{$key}[0]\t$hash{$key}[1]\t$hash{$key}[2]\t$hash{$key}[3]\t$hash{$key}[4]\n";
      }
    }
    close OUT;'
    

    echo '---
title: "scATAC pipeline, ${TestID}"
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
    echo "# Alignment state of human reference" >>${Outdir}/report/report.Rmd   
    cat ${HsAln} >> ${Outdir}/report/report.Rmd

    echo "" >>${Outdir}/report/report.Rmd
    echo "# Alignment state of mouse reference (mm9)" >>${Outdir}/report/report.Rmd
    cat ${MmAln} >> ${Outdir}/report/report.Rmd

    echo "" >>${Outdir}/report/report.Rmd
    echo "# Density of barcode reads & barcode capture effciency" >>${Outdir}/report/report.Rmd

    echo '```{r echo=FALSE, fig.height = 3,  fig.width = 12, fig.align = "center"}
library("ggplot2")
library("gridExtra")
library("data.table")
## usable reads cutoff
usable_reads <- 500
## human / mouse read rate
reads_rate <- 0.8 
## barcode count cutoff
barcode_cutoff <- 1000
## ratio of reads covered on promoter
Tss_rate <- 0.1

result <- fread("${Outdir}/report/result_${TestID}.txt", header=T)
result=subset(result,result$barcode_reads > 100)

p1 <- ggplot(data = result, mapping = aes(x = log10(barcode_reads))) + 
   geom_density(color="coral",size=1)+theme_bw()+
    geom_vline(xintercept= log10(barcode_cutoff),lty=2,col="red",lwd=1.0)+
    xlab("log10 barcode raw reads")+
    ylab("Density")+
    ggtitle("Density of barcode reads")+
       theme(plot.title = element_text(hjust = 0.5))

sub_barcode=subset(result,result$barcode_reads > 0)[,c(1,2)]
order_sub_barcode=sub_barcode[order(sub_barcode$barcode_reads,decreasing = T),]
order_sub_barcode$bar_ID=seq(1:dim(order_sub_barcode)[1])

order_sub_barcode$Evaluation=as.factor(ifelse(order_sub_barcode$barcode_reads>1000,"Barcode","Contamination"))

p2 <- ggplot(data = order_sub_barcode,aes(x = log10(bar_ID),y = log10(barcode_reads),color=Evaluation)) +
     geom_point(alpha=0.7, size=1.75)+
     theme_gray() +
     theme_bw()+
     ylab("log10 Reads number") + 
     xlab("Barcode ID") +
     ggtitle("scATAC barcode capture effciency")+
     theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
     scale_color_brewer(palette="Set2",direction = -1,labels = paste(table(order_sub_barcode$Evaluation),"  ",levels(order_sub_barcode$Evaluation)))+
     scale_x_continuous(breaks=order_sub_barcode$bar_ID, labels = 10^order_sub_barcode$bar_ID)

grid.arrange(p1,p2,ncol=2)
```

# Evaluate usable reads mapped to human reference

```{r echo=FALSE, fig.height = 3,  fig.width = 12, fig.align = "center"}
sub=subset(result,result$barcode_reads>barcode_cutoff)

sub$hs_reads_rate=sub$Hs_usable_reads/(sub$Hs_usable_reads+sub$Mm_usable_reads)
sub$mm_reads_rate=sub$Mm_usable_reads/(sub$Hs_usable_reads+sub$Mm_usable_reads)
sub$map_rate_max=apply(sub[,c(7,8)],1,max)

sub$species <- as.factor(ifelse(sub$map_rate_max > reads_rate,
                                   ifelse(sub$hs_reads_rate > sub$mm_reads_rate ,"Human","Mouse"),"Doublet"))
sub=na.omit(sub)

sub_HS=subset(sub,sub$species=="Human")

p3<- ggplot(data = sub_HS, mapping = aes(x = log10(Hs_usable_reads))) + 
  geom_density(color="slate blue",size=1)+theme_bw()+
  geom_vline(xintercept= log10(usable_reads),lty=2,col="red",lwd=0.5)+
  xlab("log10 Human usable reads")+
  ylab("Density")+
  ggtitle("Density of Hs usable reads")+
  theme(plot.title = element_text(hjust = 0.5))

p4<- ggplot(data = sub_HS, mapping = aes(x = sub_HS$Hs_reads_in_promoters/sub_HS$Hs_usable_reads)) + 
  geom_density(color="medium sea green",size=1)+theme_bw()+
  geom_vline(xintercept=Tss_rate,lty=2,col="red",lwd=0.5)+
  xlab("rate of reads in promoters")+
  ylab("Density")+
  ggtitle("Density of Human reads in promoters(4Kb)")+
  theme(plot.title = element_text(hjust = 0.5))

sub_HS$Tags=as.factor(ifelse(sub_HS$Hs_usable_reads>usable_reads & sub_HS$Hs_reads_in_promoters/sub_HS$Hs_usable_reads >Tss_rate,"Usable Hs cells","FALSE"))

p5 <- ggplot(data=sub_HS,aes(x=log10(sub_HS$Hs_usable_reads), y =sub_HS$Hs_reads_in_promoters/sub_HS$Hs_usable_reads,
                     color=Tags))+
  geom_point(alpha=0.7, size=1.75) +
  xlab("log10 Human usable reads") + 
  ylab("Human rate of reads in promoter") +
  ggtitle("Human QC")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_vline(xintercept= log10(usable_reads),lty=2,col="red",lwd=0.5)+
  geom_hline(yintercept= Tss_rate,lty=2,col="red",lwd=0.5)+
  scale_color_manual(values = c("gray","navy"),labels = paste(table(sub_HS$Tags),"  ",levels(sub_HS$Tags)))


grid.arrange(p3,p4,p5,ncol=3,widths=c(3,3,6))
```

# Evaluate usable reads mapped to mouse reference

```{r echo=FALSE, fig.height = 3,  fig.width = 12, fig.align = "center"}

sub_MM=subset(sub,sub$species=="Mouse")

p6<- ggplot(data = sub_MM, mapping = aes(x = log10(Mm_usable_reads))) + 
  geom_density(color="slate blue",size=1)+theme_bw()+
  geom_vline(xintercept= log10(usable_reads),lty=2,col="red",lwd=0.5)+
  xlab("log10 Mouse usable reads")+
  ylab("Density")+
  ggtitle("Density of Mouse usable reads")+
  theme(plot.title = element_text(hjust = 0.5))

p7<- ggplot(data = sub_MM, mapping = aes(x = sub_MM$Mm_reads_in_promoters/sub_MM$Mm_usable_reads)) + 
  geom_density(color="medium sea green",size=1)+theme_bw()+
  geom_vline(xintercept= Tss_rate,lty=2,col="red",lwd=0.5)+
  xlab("rate of reads in promoters")+
  ylab("Density")+
  ggtitle("Density of Mouse reads in promoters(4Kb)")+
  theme(plot.title = element_text(hjust = 0.5))

sub_MM$Tags=as.factor(ifelse(sub_MM$Mm_usable_reads>usable_reads & sub_MM$Mm_reads_in_promoters/sub_MM$Mm_usable_reads >Tss_rate,"Usable Mm cells","FALSE"))

p8 <- ggplot(data=sub_MM,aes(x=log10(sub_MM$Mm_usable_reads), y =sub_MM$Mm_reads_in_promoters/sub_MM$Mm_usable_reads,
                             color=Tags))+
  geom_point(alpha=0.7, size=1.75) +
  xlab("log10 Mouse usable reads") + 
  ylab("Mouse rate of reads in promoter") +
  ggtitle("Mouse QC")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_vline(xintercept= log10(as.numeric(usable_reads)),lty=2,col="red",lwd=0.5)+
  geom_hline(yintercept= Tss_rate,lty=2,col="red",lwd=0.5)+
  scale_color_manual(values = c("gray","sea green"),labels = paste(table(sub_MM$Tags),"  ",levels(sub_MM$Tags)))

grid.arrange(p6,p7,p8,ncol=3,widths=c(3,3,6))

```
' >> ${Outdir}/report/report.Rmd

echo '
# Accessibility of enhancer

![](${accEnh_hs_plot})
![](${accEnh_mm_plot})
' >> ${Outdir}/report/report.Rmd

echo '
# Fragement size distribution
    
```{r echo=FALSE, fig.height = 3,  fig.width = 12, fig.align = "center"}
fragHs <- fread("${fragHs}", header=F)
fragMm <- fread("${fragMm}", header=F)
f1 <- ggplot(data=subset(fragHs,fragHs$V1<1000), aes(V1, V2)) + geom_point() + xlab("Fragement size") + ylab("Counts")
f2 <- ggplot(data=subset(fragMm,fragMm$V1<1000), aes(V1, V2)) + geom_point() + xlab("Fragement size") + ylab("Counts")
grid.arrange(f1,f2,ncol=2)
```    
' >> ${Outdir}/report/report.Rmd
echo '    
# Summary
    
```{r echo=FALSE, fig.height = 3,  fig.width = 12, fig.align = "center"}
sub$Tags=as.factor(ifelse(sub$species=="Human",
                    ifelse(sub$Hs_usable_reads>usable_reads& sub$Hs_reads_in_promoters/sub$Hs_usable_reads>Tss_rate,"Human_OK","Human_FALSE"),
                    ifelse(sub$species=="Mouse",
                           ifelse(sub$Mm_usable_reads>usable_reads & sub$Mm_reads_in_promoters/sub$Mm_usable_reads>Tss_rate,"Mouse_OK","Mouse_FALSE"),"Doublet")))

readmax=max(sub$Hs_usable_reads,sub$Mm_usable_reads)
p9<-ggplot(data=sub,aes(x=sub$Hs_usable_reads, y =sub$Mm_usable_reads,
                     color=Tags))+
  geom_point(alpha=0.7, size=1.75) +
  xlab("Human usable fragments") + 
  ylab("Mouse usable fragments") +
  ggtitle("single cell ATAC QC")+
  theme_gray() +
  theme_bw()+
  xlim(0,readmax)+
  ylim(0,readmax)+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_brewer(palette="Set2",labels = paste(table(sub$Tags),"  ",levels(sub$Tags)))


p10<-ggplot(data=sub,aes(x=log10(sub$Hs_usable_reads), y =log10(sub$Mm_usable_reads),
                        color=Tags))+
  geom_point(alpha=0.7, size=1.75) +
  xlab("log10 Human usable fragments") + 
  ylab("log10 Mouse usable fragments") +
  ggtitle("single cell ATAC QC")+
  theme_gray() +
  theme_bw()+
  xlim(0,log10(readmax))+
  ylim(0,log10(readmax))+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_brewer(palette="Set2",labels = paste(table(sub$Tags),"  ",levels(sub$Tags)))

grid.arrange(p9,p10,ncol=2)

```
' >> ${Outdir}/report/report.Rmd
    export RSTUDIO_PANDOC=${Root}/third_party/pandoc
    ${Root}/third_party/R -e 'rmarkdown::render("${Outdir}/report/report.Rmd", output_format="html_document",output_file ="${Outdir}/report/report.html")'

>>>
  
  output {
    String html_report="${Outdir}/report/report.html"
  }
}
