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






#################################################
# ATAC-seq and CUT&Tag-seq DBR finding workflow #
#################################################
library("edgeR")
library("NOISeq")

# do TMM normalization
counttable=read.table("edgeR_input_table.txt",header=F)
mobData_data=counttable[,-1]	
tmm_table=tmm(mobData_data, long = 1000, lc = 0, k = 0)
tmm_table =round(tmm_table,digits=0)
new_TMM_count=data.frame(mobData_data[,1],tmm_table)

# prepare for edgeR
group_edgeR=factor(c(rep(0,as.numeric(args[8])),rep(1,as.numeric(args[9]))),label=c(args[6],args[7]))
# create DGE List
wtest_raw=DGEList(counts=new_TMM_count[,-1], group=group_edgeR, genes=new_TMM_count[,1])
# Filter the peaks with CPM > $cpm_value
keep <- rowSums(cpm(wtest_raw)>as.numeric(args[10])) >=sum(as.numeric(args[8])+as.numeric(args[9]))  #delete the the peaks with 0 count in all samples
wtest <- wtest_raw[keep,]
# Normalizing the data (default method is TMM)
# wtest <- calcNormFactors(wtest) 
# Estimate dispersion
wtest<-estimateCommonDisp(wtest, verbose=TRUE)
wtest<-estimateTagwiseDisp(wtest)
# calculate the CPM
wcpm=cpm(wtest$pseudo.counts)
# Differential Analysis (from the tutorial of edgeR)
wet=exactTest(wtest)
wntop=topTags(wet,n=dim(wtest)[1])
wstop=wntop$table[order(rownames(wntop$table)),]
head(wstop)
wsexp=wcpm[order(rownames(wcpm)),]
wnftab=cbind(wstop,wsexp)
wnnftab=wnftab[order(wnftab[,4]),]
wftab=wnnftab
colnames(wftab)=c("genes","logFC","logCPM","PValue","FDR",rep(args[6], as.numeric(args[8])), rep(args[7], as.numeric(args[9])))
# select peaks with log(Fold-Change) > $log2FC and FDR < $qvalue
wftab_selected=wftab[abs(wftab[,2])>as.numeric(args[12]) & wftab[,5]<as.numeric(args[11]),]
#write out the data
write.table(wftab, file = 'edgeR_DBR_full_list.txt', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(wftab_selected, file = 'edgeR_DBR_significant_only.txt', sep = '\t', quote = FALSE, row.names = FALSE)
