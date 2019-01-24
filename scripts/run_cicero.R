### Get the parameters
parser = argparse::ArgumentParser(description="Run cicero")
parser$add_argument('-P','--peak', help='the peak file')
parser$add_argument('-C','--count', help='peak raw count file')
parser$add_argument('-O','--out', help='outpath')
parser$add_argument('-S','--sample', help='sample ID')
args = parser$parse_args()

###
library("cicero")
library("reshape2")

peak=read.table(args$peak)
count=read.table(args$count,header = T)
count$peak_ID=row.names(count)
peak=peak[,c(1:4)]
peak$peak_pos=paste(peak$V1,peak$V2,peak$V3,sep = "_")
peak=peak[,c(4,5)]
colnames(peak)=c("peak_ID","peak_pos")

merge=merge(peak,count,by="peak_ID")
merge=merge[,-1]
melt=melt(merge,id.vars ="peak_pos")
colnames(melt)=c("Peak","Cell","Count")

cicero_data=melt
input_cds <- make_atac_cds(cicero_data, binarize = TRUE)
agg_cds <- aggregate_nearby_peaks(input_cds)
agg_cds <- detectGenes(agg_cds)
agg_cds <- estimateSizeFactors(agg_cds)
agg_cds <- estimateDispersions(agg_cds)

agg_cds <- reduceDimension(agg_cds,
                           max_components = 2,
                           norm_method = 'log',
                           num_dim = 6,
                           reduction_method = 'tSNE',
                           verbose = T)

agg_cds <- clusterCells(agg_cds, verbose = F)
pdf(paste0(args$out,"/","Cicero_cluster_",args$sample,".pdf"))
plot_cell_clusters(agg_cds, color_by = 'as.factor(Cluster)')
dev.off()









