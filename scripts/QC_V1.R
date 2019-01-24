### Get the parameters
parser = argparse::ArgumentParser(description="Script to plot fragment size")
parser$add_argument('-C','--barcode', help='the reads of every barcode captured')
parser$add_argument('-H_r','--Hs_use', help='human usable')
parser$add_argument('-H_t','--Hs_tss', help='human usable in promoters')
parser$add_argument('-M_r','--Mm_use', help='mouse usable')
parser$add_argument('-M_t','--Mm_tss', help='mouse usable in promoters')
parser$add_argument('-E','--bar', help='barcode reads')
parser$add_argument('-R','--rate', help='human usable reads/mouse usable reads')
parser$add_argument('-U','--use', help='number of usable reads')
parser$add_argument('-T','--tss', help='rate of usable reads in promoters')
parser$add_argument('-O','--out', help='outpath')
parser$add_argument('-B','--batch', help='batch ID')
args = parser$parse_args()
###

library("ggplot2")
library("gridExtra")
library("data.table")

#-----------------------------------------------------------------------#
#                         Hs statistics                                 #
#-----------------------------------------------------------------------#
ID_hs=fread(args$Hs_use,header = F)
colnames(ID_hs)=c("read_ID_hs","barcode_ID_hs")
Usable_hs=as.data.frame(table(ID_hs$barcode_ID_hs))
colnames(Usable_hs)=c("barcode_ID","Hs_usable_reads")
TSS_hs=fread(args$Hs_tss,header = F)
colnames(TSS_hs)=c("read_ID_hs")
merge_hs=merge(TSS_hs,unique(ID_hs),by="read_ID_hs")
pro_hs=as.data.frame(table(merge_hs$barcode_ID_hs))
colnames(pro_hs)=c("barcode_ID","Hs_reads_in_promoters")
merge_1=merge(Usable_hs,pro_hs,by="barcode_ID",all = T)
#-----------------------------------------------------------------------#
#                         Mm statistics                                 #
#-----------------------------------------------------------------------#
ID_mm=fread(args$Mm_use,header = F)
colnames(ID_mm)=c("read_ID_mm","barcode_ID_mm")
Usable_mm=as.data.frame(table(ID_mm$barcode_ID_mm))
colnames(Usable_mm)=c("barcode_ID","Mm_usable_reads")
TSS_mm=fread(args$Mm_tss,header = F)
colnames(TSS_mm)=c("read_ID_mm")
merge_mm=merge(TSS_mm,unique(ID_mm),by="read_ID_mm")
pro_mm=as.data.frame(table(merge_mm$barcode_ID_mm))
colnames(pro_mm)=c("barcode_ID","Mm_reads_in_promoters")
merge_2=merge(Usable_mm,pro_mm,by="barcode_ID",all = T)
#-----------------------------------------------------------------------#
#                         Merge result                                  #
#-----------------------------------------------------------------------#
all=fread(args$barcode,header = F)
colnames(all)=c("barcode_ID","barcode_reads")
all$barcode_reads=all$barcode_reads*2
res_1=merge(all,merge_1,by="barcode_ID",all=T)
res_2=merge(res_1,merge_2,by="barcode_ID",all=T)
res_2[is.na(res_2)]=0
result=res_2
#result=subset(result,result$barcode_reads > 100)
write.table(result,paste0(args$out,"/","result-",args$batch,".txt"),sep = "\t",row.names = FALSE,quote = FALSE)
#=======================================================================#
#                              Plot                                     #
#=======================================================================#

##################plot-1 Evaluate barcode capture effciency###############
result=subset(result,result$barcode_reads > 100)

p1<- ggplot(data = result, mapping = aes(x = log10(barcode_reads))) + 
  geom_density(color="coral",size=1)+theme_bw()+
  geom_vline(xintercept= log10(as.numeric(args$bar)),lty=2,col="red",lwd=1.0)+
  xlab("log10 barcode raw reads")+
  ylab("Density")+
  ggtitle("Density of barcode reads")+
  theme(plot.title = element_text(hjust = 0.5))
#print(p1)

sub_barcode=subset(result,result$barcode_reads > 0)[,c(1,2)]
order_sub_barcode=sub_barcode[order(sub_barcode$barcode_reads,decreasing = T),]
order_sub_barcode$bar_ID=seq(1:dim(order_sub_barcode)[1])

order_sub_barcode$Evaluation=as.factor(ifelse(order_sub_barcode$barcode_reads>as.numeric(args$bar),"Barcode","Contamination"))

p2<- ggplot(data = order_sub_barcode,aes(x = log10(bar_ID),y = log10(barcode_reads),color=Evaluation)) +
     geom_point(alpha=0.7, size=1.75)+
     theme_gray() +
     theme_bw()+
     ylab("log10 Reads number") + 
     xlab("Barcode ID") +
     ggtitle("Drop-ATAC barcode carpure effciency")+
     theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
     scale_color_brewer(palette="Set2",direction = -1,labels = paste(table(order_sub_barcode$Evaluation),"  ",levels(order_sub_barcode$Evaluation)))+
     scale_x_continuous(breaks=order_sub_barcode$bar_ID, labels = 10^order_sub_barcode$bar_ID)
#print(p2)

pdf(paste0(args$out,"/","plot-1-Barcode-",args$batch,".pdf"),width = 13,height = 3.5)
grid.arrange(p1,p2,ncol=2)
dev.off()

##################plot-2 Evaluate Hs usable cell ###################

sub=subset(result,result$barcode_reads>as.numeric(args$bar))

ma=as.numeric(args$rate)

sub$hs_reads_rate=sub$Hs_usable_reads/(sub$Hs_usable_reads+sub$Mm_usable_reads)
sub$mm_reads_rate=sub$Mm_usable_reads/(sub$Hs_usable_reads+sub$Mm_usable_reads)
sub$map_rate_max=apply(sub[,c(7,8)],1,max)

sub$species <- as.factor(ifelse(sub$map_rate_max > ma,
                                   ifelse(sub$hs_reads_rate > sub$mm_reads_rate ,'Human','Mouse'),'Doublet'))
sub=na.omit(sub)

sub_HS=subset(sub,sub$species=="Human")

p3<- ggplot(data = sub_HS, mapping = aes(x = log10(Hs_usable_reads))) + 
  geom_density(color="slate blue",size=1)+theme_bw()+
  geom_vline(xintercept= log10(as.numeric(args$use)),lty=2,col="red",lwd=0.5)+
  xlab("log10 Human usable reads")+
  ylab("Density")+
  ggtitle("Density of Hs usable reads")+
  theme(plot.title = element_text(hjust = 0.5))
#print(p3)

p4<- ggplot(data = sub_HS, mapping = aes(x = sub_HS$Hs_reads_in_promoters/sub_HS$Hs_usable_reads)) + 
  geom_density(color="medium sea green",size=1)+theme_bw()+
  geom_vline(xintercept= as.numeric(args$tss),lty=2,col="red",lwd=0.5)+
  xlab("rate of reads in promoters")+
  ylab("Density")+
  ggtitle("Density of Human reads in promoters(4Kb)")+
  theme(plot.title = element_text(hjust = 0.5))
#print(p4)

sub_HS$Tags=as.factor(ifelse(sub_HS$Hs_usable_reads>as.numeric(args$use) & sub_HS$Hs_reads_in_promoters/sub_HS$Hs_usable_reads >as.numeric(args$tss),"Usable Hs cells","FALSE"))

p5 <- ggplot(data=sub_HS,aes(x=log10(sub_HS$Hs_usable_reads), y =sub_HS$Hs_reads_in_promoters/sub_HS$Hs_usable_reads,
                     color=Tags))+
  geom_point(alpha=0.7, size=1.75) +
  xlab("log10 Human usable reads") + 
  ylab("Human rate of reads in promoter") +
  ggtitle("Human QC")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_vline(xintercept= log10(as.numeric(args$use)),lty=2,col="red",lwd=0.5)+
  geom_hline(yintercept= as.numeric(args$tss),lty=2,col="red",lwd=0.5)+
  scale_color_manual(values = c("gray","navy"),labels = paste(table(sub_HS$Tags),"  ",levels(sub_HS$Tags)))

pdf(paste0(args$out,"/","plot-2-Evaluate-HS-",args$batch,".pdf"),width = 16,height = 3.5)
grid.arrange(p3,p4,p5,ncol=3)
dev.off()

##################plot-3 Evaluate Mm usable cell ###################
sub_MM=subset(sub,sub$species=="Mouse")

p6<- ggplot(data = sub_MM, mapping = aes(x = log10(Mm_usable_reads))) + 
  geom_density(color="slate blue",size=1)+theme_bw()+
  geom_vline(xintercept= log10(as.numeric(args$use)),lty=2,col="red",lwd=0.5)+
  xlab("log10 Mouse usable reads")+
  ylab("Density")+
  ggtitle("Density of Mouse usable reads")+
  theme(plot.title = element_text(hjust = 0.5))
#print(p6)

p7<- ggplot(data = sub_MM, mapping = aes(x = sub_MM$Mm_reads_in_promoters/sub_MM$Mm_usable_reads)) + 
  geom_density(color="medium sea green",size=1)+theme_bw()+
  geom_vline(xintercept= as.numeric(args$tss),lty=2,col="red",lwd=0.5)+
  xlab("rate of reads in promoters")+
  ylab("Density")+
  ggtitle("Density of Mouse reads in promoters(4Kb)")+
  theme(plot.title = element_text(hjust = 0.5))
#print(p7)

sub_MM$Tags=as.factor(ifelse(sub_MM$Mm_usable_reads>as.numeric(args$use) & sub_MM$Mm_reads_in_promoters/sub_MM$Mm_usable_reads >as.numeric(args$tss),"Usable Mm cells","FALSE"))

p8 <- ggplot(data=sub_MM,aes(x=log10(sub_MM$Mm_usable_reads), y =sub_MM$Mm_reads_in_promoters/sub_MM$Mm_usable_reads,
                             color=Tags))+
  geom_point(alpha=0.7, size=1.75) +
  xlab("log10 Mouse usable reads") + 
  ylab("Mouse rate of reads in promoter") +
  ggtitle("Mouse QC")+
  theme_gray() +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  geom_vline(xintercept= log10(as.numeric(args$use)),lty=2,col="red",lwd=0.5)+
  geom_hline(yintercept= as.numeric(args$tss),lty=2,col="red",lwd=0.5)+
  scale_color_manual(values = c("gray","sea green"),labels = paste(table(sub_MM$Tags),"  ",levels(sub_MM$Tags)))
#print(p8)

pdf(paste0(args$out,"/","plot-3-Evaluate-MM-",args$batch,".pdf"),width = 16,height = 3.5)
grid.arrange(p6,p7,p8,ncol=3)
dev.off()

##################plot-4 Summary ###################
sub$Tags=as.factor(ifelse(sub$species=="Human",
                    ifelse(sub$Hs_usable_reads> as.numeric(args$use)& sub$Hs_reads_in_promoters/sub$Hs_usable_reads>as.numeric(args$tss),"Human_OK","Human_FALSE"),
                    ifelse(sub$species=="Mouse",
                           ifelse(sub$Mm_usable_reads>as.numeric(args$use) & sub$Mm_reads_in_promoters/sub$Mm_usable_reads>as.numeric(args$tss),"Mouse_OK","Mouse_FALSE"),"Doublet")))

readmax=max(sub$Hs_usable_reads,sub$Mm_usable_reads)
p9<-ggplot(data=sub,aes(x=sub$Hs_usable_reads, y =sub$Mm_usable_reads,
                     color=Tags))+
  geom_point(alpha=0.7, size=1.75) +
  xlab("Human usable fragments") + 
  ylab("Mouse usable fragments") +
  ggtitle("Droplet ATAC QC")+
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
  ggtitle("Droplet ATAC QC")+
  theme_gray() +
  theme_bw()+
  xlim(0,log10(readmax))+
  ylim(0,log10(readmax))+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_color_brewer(palette="Set2",labels = paste(table(sub$Tags),"  ",levels(sub$Tags)))
#print(p10)

pdf(paste0(args$out,"/","plot-4-Summary-",args$batch,".pdf"),width = 12,height = 4.5)
grid.arrange(p9,p10,ncol=2)
dev.off()








