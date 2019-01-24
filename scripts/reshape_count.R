# author: Xiaoyu Wei

### Get the parameters
parser = argparse::ArgumentParser(description="Run cicero")
parser$add_argument('-I','--input', help='the count file')
parser$add_argument('-o','--out', help='output file')
parser$add_argument('-S','--sample', help='sample ID')
args = parser$parse_args()

###
library(data.table)

input=fread(args$input)

input$peak=paste(input$V1,input$V2,input$V3,sep = "_")
input=input[,c(6,5)]
colnames(input)=c("peak","barcode")

input$V=1

out <- input[,list(count=sum(V)),by=list(peak,barcode)]
out1 <- reshape(out,direction='wide',idvar='peak', timevar='barcode')
out2 <- as.data.frame(out1)
out2[is.na(out2)]=0

names(out2) <- gsub("count.", "", names(out2))
write.table(out2,args$out,sep="\t",quote=FALSE,row.names=FALSE)
