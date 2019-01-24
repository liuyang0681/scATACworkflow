#### args 

args = commandArgs(TRUE)

infile = args[1]
#bed = args[2]
                                        #out.prefix = args[3]
out.prefix = args[2]

#### library
library(chromVAR)
library(GenomicRanges)
library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm9)
library(motifmatchr)
library(TFBSTools)
library(MotifDb)
library(JASPAR2016)
library(universalmotif)
library(data.table)
library(ggplot2)
library(pheatmap)
library(viridis)
library(ggrepel)
library(gridExtra)
library(densityClust)
library(BiocParallel)

#### Function
plotTFvariation <- function(variability.atac.counts.filter,n = 20){
  res_df <- cbind(variability.atac.counts.filter, rank = rank(-1 * variability.atac.counts.filter$variability, ties.method = "random"), 
                  annotation = variability.atac.counts.filter$name)
  
  top_df <- res_df[res_df$rank <= n, ]
  
  p <- ggplot(res_df, aes_string(x = "rank", y = "variability", min = "bootstrap_lower_bound",
                                 max = "bootstrap_upper_bound", label = "annotation")) + 
    geom_errorbar(colour = "grey") + 
    geom_point(colour = "purple",size = 1) + 
    xlab("Sorted TFs") + 
    ylab("Variability") + 
    theme_bw() + 
    scale_y_continuous(expand= c(0, 0),limits = c(0, max(res_df$bootstrap_upper_bound, na.rm = TRUE) * 1.05)) + 
    geom_text_repel(data = top_df, size = 3, col = "Black") # require ggrepel package
  return(p)
}

#### set core number
register(MulticoreParam(24))

##### Data analysis and Figures

#source("../Desktop/plotTFvariation.R")
#count = fread(args[1])
#promoter = fread(args[2])

print("Read RDS")

data = read.table(infile,header=T,row.names=1,sep="\t")

Sum_peak = rowSums(data)
Sum_ReadNum = colSums(data)

data = data[Sum_peak >=10,Sum_ReadNum>=1000]

print("SummariExpressiment")
peak.bed = do.call(rbind,strsplit(x=rownames(data),split="_"))
peak.bed = as.data.frame(peak.bed) 
#peak.bed = read.table(bed, sep="\t",col.names=c("chr","start","end","name","var1","var2","var3","var4"))
#peak.bed = peak.bed[peak.bed$name %in% rownames(data),]
colnames(peak.bed) = c("chr","start","end")
peak.gr <- makeGRangesFromDataFrame(peak.bed, keep.extra.columns = T)
#peak = paste0(peak.bed[,1],"_",peak.bed[,2],"_",peak.bed[,3])
names(peak.gr) = rownames(data)


# 
data = as.data.frame(data)

sample.info = as.data.frame(cbind(Sample=colnames(data), Batch =sub("_\\d+$","",colnames(data))))

atac.macaca <- SummarizedExperiment(assays = list(counts = as.matrix(data)),
                                rowRanges = peak.gr,colData = sample.info)

mf = BSgenome.Mmusculus.UCSC.mm9

seqnames(mf) = sub("MFA","chr",seqnames(mf))

atac.macaca <- addGCBias(atac.macaca, genome = mf)

atac.macaca.filter <- filterPeaks(atac.macaca, non_overlapping = F)

motifs <- getJasparMotifs()
#matrices.TF <- query(MotifDb, c('hsapiens'))
#matrices.TF <- query(MotifDb, c('hsapiens', 'mmusculus'))
#matrices.TF_v2 <- convert_motifs(matrices.TF, "TFBSTools-PFMatrix")
#PFMatrixList <- convert_motifs(motifs, "TFBSTools-PFMatrix")
#motifs <- do.call(PFMatrixList, matrices.TF_v2)

atac.macaca.filter.motif_ix_v2 <- matchMotifs(motifs, atac.macaca.filter, 
                                              genome = mf)

#atac.macaca.filter.motif_ix_v2_2 = atac.macaca.filter.motif_ix_v2[,-c(grep("FOS|JUN|JDP",
#                                                                           colnames(atac.macaca.filter.motif_ix_v2)))]

dev.atac.macaca.filter <- computeDeviations(object = atac.macaca.filter, 
                                            annotations = atac.macaca.filter.motif_ix_v2)

variability.atac.macaca.filter <- computeVariability(dev.atac.macaca.filter)
#variability.atac.macaca.filter2 = variability.atac.macaca.filter[-c(grep("FOS|JUN|JDP",variability.atac.macaca.filter$name)),]

pdf(paste0(out.prefix,".TF_deviation.pdf"))
tfv=plotTFvariation(variability.atac.macaca.filter)
print(tfv)

set.seed(1)
tsne_results <- deviationsTsne(dev.atac.macaca.filter, threshold = 1.03, perplexity = 50)

#tsne_plots <- plotDeviationsTsne(dev.atac.macaca.filter, tsne_results, annotation_name = "LMX1B", shiny = FALSE)

TFs = as.vector(variability.atac.macaca.filter[order(variability.atac.macaca.filter$variability,decreasing=T),]$name[1:5])

for(i in 1:length(TFs)){
  pdf(paste0(out.prefix,".TF_",i,".pdf"))
  name = TFs[i]
  a_plots <- plotDeviationsTsne(dev.atac.macaca.filter, tsne_results, 
                                annotation_name = name, shiny = FALSE)
  print(a_plots)
}
dev.off()


pdf(paste0(out.prefix,".Cell_Cluster.pdf"))

dis <- dist(tsne_results)
denCluster = densityClust(dis, gaussian=TRUE)

#plot(denCluster)

rho_cutoff <- median(denCluster$rho)
delta_cutoff <- round(denCluster$dc + 0.5)

denClust <- findClusters(denCluster, 
                          rho= rho_cutoff, 
                          delta= delta_cutoff)

clusters <- denClust$cluster
plot(tsne_results, col=clusters, cex=0.5, pch=19)

#a_plots <- plotDeviationsTsne(dev.atac.macaca.filter, tsne_results,
#                              sample_column = "Batch", shiny = FALSE)
print(a_plots)

dev.off()

save.image(paste0(out.prefix,".TF_chromVAR.Rdata"))
