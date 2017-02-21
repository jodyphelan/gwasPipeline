library(scales)

args = commandArgs(trailingOnly=TRUE)

pca<-read.table(args[1])
vars<-read.table(args[2])
meta<-read.table(args[3])

meta<-meta[match(pca$V2,meta$V1),]

meta.uniq<-unique(meta$V2)
cols.uniq<-alpha(rainbow(length(meta.uniq)),0.5)
cols<-cols.uniq[match(meta$V2,meta.uniq)]

xlabs<-paste("PC1 (",vars$V1[1],"%)",sep="")
ylabs<-paste("PC1 (",vars$V1[2],"%)",sep="")
png(args[4],800,600)
plot(pca$V3,pca$V4,col=cols,pch=20,cex=2,xlab=xlabs,ylab=ylabs)
legend("bottomleft",legend=meta.uniq,fill=cols.uniq)
dev.off()
