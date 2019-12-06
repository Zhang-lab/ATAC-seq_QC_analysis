#!/usr/bin/env Rscript
args=commandArgs()

#################################
# ATAC-seq DOR finding workflow #
#################################
library("DESeq2")

### load count matrix
countdata=read.table("full.count",header=F)
rownames(countdata)=paste0(rep("merged_peak",nrow(countdata)),seq(1,nrow(countdata)))

group=factor(c(rep(0,as.numeric(args[8])),rep(1,as.numeric(args[9]))),label=c(args[6],args[7]))
colData=data.frame(group)
rownames(colData)=colnames(countdata)[-(1:3)]

dds=DESeqDataSetFromMatrix(countData=countdata[,-(1:3)],colData=colData,design=~group)

dds=dds[rowSums(fpm(dds,robust=F)>1)>=2,]

# differential open analysis
dds=DESeq(dds)
res=results(dds,contrast=c("group",args[6],args[7]),alpha=0.05,independentFiltering=F,cooksCutoff=Inf)

open=cbind(countdata[rownames(countdata)%in%rownames(assay(dds)),1:3],rownames(assay(dds)),res$log2FoldChange,res$padj)
colnames(open)=c("chr","start","end","peak_index","log2FoldChange","padj")
write.table(open,file="DOR_all_regions.bed",sep="\t",row.names=F,col.names=T,quote=FALSE)
open=open[abs(open[,5])>1&open[,6]<0.01,]
write.table(open,file="DOR_full.bed",sep="\t",row.names=F,col.names=T,quote=FALSE)
