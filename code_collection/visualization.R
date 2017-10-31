#!/usr/bin/env Rscript
args=commandArgs()

name=args[6]
name=paste(name,"result",sep=".")
refpath=args[7]

test=try(library(ggplot2),silent=T)
if(class(test)=='try-error') {
	install.packages('ggplot2')
	library(ggplot2)
}
test=try(library(cowplot),silent=T)
if(class(test)=='try-error') {
	install.packages('cowplot')
	library(cowplot)
}

# chromosome distribution
chr=read.table(paste("chrom_count",name,sep="_"))
r=as.numeric(rownames(chr[which(chr$V1=='chrM'),]))
chr=rbind(chr[-r,],chr[r,])
rownames(chr)=NULL
chr$V1=factor(c(nrow(chr):1),labels=rev(as.character(chr$V1)))
chr$V2=chr$V2/sum(chr$V2)
chr$V3=chr$V3/sum(chr$V3)
colnames(chr)=c("chromosome","Total mapped reads","Non-redundant uniquely mapped reads")
chr=data.frame(chr[1],stack(chr[2:3]))

png("plot2.2_reads_distri_in_chrom.png",height=640,width=6600,res=300)
ggplot(chr,aes(x=ind,y=values,fill=chromosome))+
	ggtitle("Stacked barplotm of reads percentage in each chromosome")+
	geom_bar(stat="identity",color="black")+
	scale_y_continuous(name="Percentage of reads")+
	scale_x_discrete(name="")+
	guides(fill=guide_legend(nrow=1,byrow=TRUE,reverse=T))+
	theme_bw()+theme_classic()+coord_flip()+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"),
		legend.text=element_text(size=10,face="bold"),
		legend.title=element_text(size=10,face="bold"),
		legend.position="bottom")
dev.off()

# saturation analysis
saturate=read.table(paste("saturation",name,sep="_"),header=T,sep='\t')
saturate=saturate[,-2]
colnames(saturate)=c("depth","peak","percentage","marker")

png("plot4.5_saturation.png",height=2300,width=2800,res=300)
plot=ggplot(saturate,aes(x=depth,y=100*percentage))+
	geom_point()+geom_line(size=1)+expand_limits(x=0,y=0)+
	geom_text(aes(label=peak),hjust=0,vjust=1.4,family="Tahoma",fontface=2,size=3)+
	ggtitle("Line chart of recaptured coverage ratio labeled with \n numbers of peaks by subsampling original library")+
	scale_x_continuous(name="Percentage of original library size",breaks=seq(0,100,by=10))+
	geom_vline(xintercept=50,linetype="dotted",size=1)+
	theme_bw()+theme_classic()+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"))
if(max(saturate$percentage)<=1) {
	plot+scale_y_continuous(name="Percentage of coverage ratio in original peaks",limits=c(0,100))
} else {
	plot+scale_y_continuous(name="Percentage of coverage ratio in original peaks")
}
dev.off()

# peak legth distribution
peakl=read.table(paste("peak_length_distri",name,sep="_"))
peakl=rbind(peakl[which(peakl$V1<1500),],c(1500,sum(peakl[which(peakl$V1>=1500),2])))
peakl$V2=peakl$V2/sum(peakl$V2)

png("plot3.3_peak_length.png",height=2000,width=3000,res=300)
ggplot(peakl,aes(x=V1,y=..scaled..,weight=V2))+geom_density(size=1,adjust=0.2)+
	ggtitle("Density plot of peaks length distribution (Adjust=0.2)")+
	theme_bw()+theme_classic()+expand_limits(x=0,y=0)+
	scale_y_continuous(name="Density",limits=c(0,1))+
	scale_x_continuous(name="Length of peaks",breaks=seq(0,max(peakl$V1),by=150))+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"))
dev.off()

# insertion distribution
insert=read.table(paste("insertion_distri",name,sep="_"))
insert$V2=insert$V2/sum(insert$V2)

png("plot3.1_insertion_size.png",height=1800,width=2400,res=300)
ggplot(insert,aes(x=V1,y=..scaled..,weight=V2))+geom_density(size=1,adjust=0.2)+
	ggtitle("Density plot of insertion size distribution (Adjust=0.2)")+
	theme_bw()+theme_classic()+
	scale_y_continuous(name="Density",limits=c(0,1))+
	scale_x_continuous(name="Insertion size",limits=c(0,500))+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"))
dev.off()

# dedup
ref.dedup=read.table(paste(refpath,'merged_dedup_percentage.txt',sep='/'),sep='\t',header=T) ####
i=1
filename=c()
ratio=c() 
while(i<=dim(ref.dedup)[1]) {
	filename=c(filename,strsplit(as.character(ref.dedup[i,1]),'_1',fixed=T)[[1]][1])
	ratio=c(ratio,100-(ref.dedup[i,2]+ref.dedup[i+1,2])/2)
	i=i+2
}
ref.dedup=data.frame(filename,ratio,class="ENCODE PE")
dedup=read.table(paste("dedup_percentage",name,sep="_"),header=T,sep='\t')
if(dim(dedup)[1]>1) {
	dedup=data.frame(filename=strsplit(as.character(dedup[1,1]),'_1',fixed=T)[[1]][1],ratio=100-(dedup[1,2]+dedup[2,2])/2,class='Sample')
	data_type="Paired-end data"
} else {
	dedup=data.frame(filename=strsplit(as.character(dedup[1,1]),'_1',fixed=T)[[1]][1],ratio=100-dedup[1,2],class='Sample')
	data_type="Single-end data"
}
ref.dedup=rbind(ref.dedup,dedup)

# enrichment
ref.enrich2=read.table(paste(refpath,'merged_coding_promoter_peak_enrichment.txt',sep='/'),header=T,sep='\t') ####
enrich2=read.table(paste("enrichment_ratio_in_promoter",name,sep="_"),header=T,sep='\t')
ref.enrich2=data.frame(enrichment_ratio=ref.enrich2[,5],class="ENCODE PE")
enrich2=data.frame(enrichment_ratio=enrich2[,5],class="Sample")
ref.enrich2=rbind(ref.enrich2,enrich2)
ref.enrich2[,3]="Enrichment ratio in coding promoter regions"

ref.enrich3=read.table(paste(refpath,'merged_new_enrichment.txt',sep='/'),header=T,sep='\t') ####
enrich3=read.table(paste("new_enrichment",name,sep="_"),header=T,sep='\t')
ref.enrich3=data.frame(enrichment_ratio=ref.enrich3[,6],class="ENCODE PE")
enrich3=data.frame(enrichment_ratio=enrich3[,5],class="Sample")
ref.enrich3=rbind(ref.enrich3,enrich3)
ref.enrich3[,3]="Normalized enrichment ratio"

test=rbind(ref.enrich2,ref.enrich3)

png("plot4.2.2_peaks_enrichment_ratio.png",height=1800,width=2800,res=300)
options(warn=-1)
plot=ggplot(test,aes(x=class,y=enrichment_ratio,fill=class))+
	stat_boxplot(geom="errorbar",size=1,width=0.3,aes(colour=class))+
    geom_boxplot(outlier.shape=NA,width=0.3,lwd=1,fatten=1,aes(colour=class))+
    scale_x_discrete(name="")+
	ggtitle("Boxplot of peaks enrichment ratio")+
	scale_fill_manual(values=c(`ENCODE PE`="grey",Sample="red"))+
	scale_colour_manual(values=c(`ENCODE PE`="black",Sample="red"))+
	facet_grid(.~V3)+theme(legend.position="none")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		strip.text.x=element_text(size=12,face="bold"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold",colour=c(`ENCODE PE`="black",Sample="red")),
		axis.text.y=element_text(size=10,face="bold"))
if(enrich3[1,1]<130) {
	plot+scale_y_continuous(name="Enrichment ratio",limits=c(0,100))
} else {
	plot+scale_y_continuous(name="Enrichment ratio")
}
options(warn=0)
dev.off()

# mapping status
ref.map=read.table(paste(refpath,'merged_mapping_status.txt',sep='/'),header=T,sep='\t') ####
ref.map=ref.map[,c(2,6,10,12)]
ref.map=data.frame(ref.map,class="ENCODE PE")
map=read.table(paste("mapping_status",name,sep="_"),header=T,sep='\t')
map=map[,c(2,6,10,12)]
map=data.frame(map,class="Sample")
ref.map=rbind(ref.map,map)
ref.map$nodup_ratio=1-ref.map$nodup_ratio
map$nodup_ratio=1-map$nodup_ratio

png("plot4.3_PCR_duplicates_percentage.png",height=2000,width=1400,res=300)
options(warn=-1)
plot=ggplot(ref.map,aes(x=class,y=100*nodup_ratio,fill=class))+
	stat_boxplot(geom="errorbar",size=1,width=0.3,aes(colour=class))+
    geom_boxplot(outlier.shape=NA,width=0.3,lwd=1,fatten=1,aes(colour=class))+
    scale_x_discrete(name="")+
	ggtitle("Boxplot of PCR duplicates percentage")+
	scale_fill_manual(values=c(`ENCODE PE`="grey",Sample="red"))+
	scale_colour_manual(values=c(`ENCODE PE`="black",Sample="red"))+
	theme_bw()+theme_classic()+theme(legend.position="none")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold",colour=c(`ENCODE PE`="black",Sample="red")),
		axis.text.y=element_text(size=10,face="bold"))
if(map$nodup_ratio*100<30) {
	plot+scale_y_continuous(name="Percentage of PCR duplicates",limits=c(0,30))
} else {
	plot+scale_y_continuous(name="Percentage of PCR duplicates")
}
options(warn=0)
dev.off()

png("plot4.1_RUP.png",height=2000,width=1400,res=300)
options(warn=-1)
plot=ggplot(ref.map,aes(x=class,y=rup_ratio,fill=class))+
	stat_boxplot(geom="errorbar",size=1,width=0.3,aes(colour=class))+
    geom_boxplot(outlier.shape=NA,width=0.3,lwd=1,fatten=1,aes(colour=class))+
    scale_x_discrete(name="")+
	ggtitle("Boxplot of reads under peaks ratio")+
	scale_fill_manual(values=c(`ENCODE PE`="grey",Sample="red"))+
	scale_colour_manual(values=c(`ENCODE PE`="black",Sample="red"))+
	theme_bw()+theme_classic()+theme(legend.position="none")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold",colour=c(`ENCODE PE`="black",Sample="red")),
		axis.text.y=element_text(size=10,face="bold"))
if(map$rup_ratio<50) {
	plot+scale_y_continuous(name="Reads under peaks ratio",limits=c(0,50))
} else {
	plot+scale_y_continuous(name="Reads under peaks ratio")
}
options(warn=0)
dev.off()

colnames(ref.map)[1:2]=c("Total reads","Non-redundant uniquely mapped reads")
test=data.frame(ref.map[5],stack(ref.map[1:2]),row.names=NULL)

png("plot3.1_library_reads_distri.png",height=1800,width=2800,res=300)
plot=ggplot(test,aes(x=class,y=values,fill=class))+
	stat_boxplot(geom="errorbar",size=1,width=0.3,aes(colour=class))+
    geom_boxplot(outlier.shape=NA,width=0.3,lwd=1,fatten=1,aes(colour=class))+
    scale_x_discrete(name="")+
	ggtitle("Boxplot of reads distribution")+
	scale_fill_manual(values=c(`ENCODE PE`="grey",Sample="red"))+
	scale_colour_manual(values=c(`ENCODE PE`="black",Sample="red"))+
	facet_grid(.~ind)+theme(legend.position="none")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		strip.text.x=element_text(size=12,face="bold"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold",colour=c(`ENCODE PE`="black",Sample="red")),
		axis.text.y=element_text(size=10,face="bold"))
if(map[1,1]<10e+7) {
	plot+scale_y_continuous(name="Number of reads",limits=c(0,10e+7))
} else {
	plot+scale_y_continuous(name="Number of reads")
}
dev.off()

# Yield plot
yield=read.table(paste("yield",name,sep="_"),sep='\t',header=T)
yield=yield[yield$TOTAL_READS<=1e8,]

png("plot2.3_yield_distinction.png",height=1800,width=2600,res=300)
ggplot(yield,aes(x=TOTAL_READS,y=EXPECTED_DISTINCT))+
	geom_ribbon(aes(ymin=yield[,3],ymax=yield[,4],x=TOTAL_READS,fill="Range of confidence interval"),alpha=0.5)+
	geom_point()+geom_line(size=1)+
	geom_point(aes(x=map[,1],y=map[,2],color="red"),show.legend=F)+
	ggtitle("Line chart of yield distinction")+
	scale_x_continuous(name="Total reads")+
	scale_y_continuous(name="Expected distinction")+
	theme_bw()+theme_classic()+
	scale_fill_manual("Reference ribbon",values="grey")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.6),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"),
		legend.text=element_text(size=8,face="bold"),
		legend.title=element_text(size=12,face="bold"),
		legend.position="bottom")
dev.off()

# promoter
promoter=read.table(paste("promoter_percentage",name,sep="_"),header=T,sep='\t')
peak=data.frame(group=c("Peaks in promoter region","Peaks in non-promoter region"),value=as.numeric(promoter[1,1:2]))
read=data.frame(group=c("Reads under peaks in promoter region","Reads under peaks in non-promoter region"),value=as.numeric(promoter[1,3:4]))

blank_theme=theme_minimal()+theme(
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	panel.border=element_blank(),
	panel.grid=element_blank(),
	axis.ticks=element_blank(),
	plot.title=element_text(size=14, face="bold"))

p1=ggplot(peak,aes(x="",y=value,fill=group))+
	geom_bar(width=1,stat="identity")+
	coord_polar("y",start=0)+
	ggtitle("Pie chart of peaks distribution")+
	scale_fill_grey()+blank_theme+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,face="bold",family="Tahoma"),
		axis.text.x=element_blank())+
	geom_text(aes(y=value/2+c(0,cumsum(value)[-length(value)]),label=value),size=5,family="Tahoma",fontface="bold")

p2=ggplot(read,aes(x="",y=value,fill=group))+
	geom_bar(width=1,stat="identity")+
	coord_polar("y",start=0)+
	ggtitle("Pie chart of reads distribution")+
	scale_fill_grey()+blank_theme+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,face="bold",family="Tahoma"),
		axis.text.x=element_blank())+
	geom_text(aes(y=value/2+c(0,cumsum(value)[-length(value)]),label=value),size=5,family="Tahoma",fontface="bold")

png("plot4.6_promoter-peak_count.png",height=1500,width=5000,res=300)
options(warn=-1)
plot_grid(p2,p1,ncol=2,nrow=1,rel_widths=c(1,0.88))
options(warn=0)
dev.off()

top=read.table(paste("bin",name,sep="_"),sep='\t')
top=data.frame(rank=as.numeric(top$V1),index=as.numeric(top$V2))

png("plot4.6_promoter_distribution_among_peaks.png",height=1800,width=3200,res=300)
ggplot(top,aes(sort(rank),index))+geom_density(stat="identity")+
	ggtitle("Percentage of promoter regions")+
	theme_bw()+theme_classic()+
	scale_y_continuous(name="Percentage of promoter regions",limits=c(0,1))+
	scale_x_continuous(name="Percentage of Top peaks")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=8,face="bold"),
		axis.text.y=element_text(size=8,face="bold"))
dev.off()

# TXT report
name=args[6]
report=list(Library=name)
name=paste(name,"result",sep=".")

genome=args[8]
report=append(report,list(Pipeline.version="V1",Genome=genome,Data.type=data_type))

ref.useful=read.table(paste(refpath,'merged_useful_reads.txt',sep='/'),header=T,sep='\t')
useful=read.table(paste("useful_reads",name,sep="_"),header=T,sep='\t')

ref.map=read.table(paste(refpath,'merged_mapping_status.txt',sep='/'),header=T,sep='\t') ####
refer=c(paste(round(mean(ref.map$total)),'(SD:',round(sd(ref.map$total)),')',sep=''))
refer=c(refer,paste(round(mean(ref.map$mapped)),'(SD:',round(sd(ref.map$mapped)),')',sep=''))
refer=c(refer,paste(round(mean(ref.map[,6])),'(SD:',round(sd(ref.map[,6])),')',sep=''))

map=read.table(paste("mapping_status",name,sep="_"),header=T,sep='\t')
samples=as.numeric(c(map$total,map$mapped,map[,6]))

ref.chr=read.table(paste(refpath,'merged_chrom_count.txt',sep='/'),header=T,sep='\t')
rownames(ref.chr)=ref.chr$chrom
ref.chr=ref.chr[,-1]
ref.chr=as.matrix(ref.chr[,seq(4,269,5)])

chr=read.table(paste("chrom_count",name,sep="_"),sep='\t')
if(genome!="danRer10") {
	refer=c(refer,paste(round(mean(ref.chr[20,]/ref.map[,6]),4)*100,'%(SD:',round(sd(ref.chr[20,]/ref.map[,6]),4)*100,'%)',sep=''))
	refer=c(refer,paste(round(mean(ref.chr[21,]/ref.map[,6]),4)*100,'%(SD:',round(sd(ref.chr[21,]/ref.map[,6]),4)*100,'%)',sep=''))
	refer=c(refer,paste(round(mean(ref.chr[22,]/ref.map[,6]),4)*100,'%(SD:',round(sd(ref.chr[22,]/ref.map[,6]),4)*100,'%)',sep=''))
	refer=c(refer,paste(round(mean(ref.map[,8])),'(SD:',round(sd(ref.map[,8])),')',sep=''))
	refer=c(refer,paste(round(mean(ref.useful[,5])),'(SD:',round(sd(ref.useful[,5])),')',sep=''))

	samples=c(samples,paste(round(chr[which(chr[,1]=='chrM'),3]/map[,6],4)*100,'%',sep=''),paste(round(chr[which(chr[,1]=='chrX'),3]/map[,6],4)*100,'%',sep=''),paste(round(chr[which(chr[,1]=='chrY'),3]/map[,6],4)*100,'%',sep=''),map[,8],useful[,5])
	library=data.frame(samples,refer)
	colnames(library)=c("Sample","ENCODE PE")
	rownames(library)=c("Total reads","Mapped reads","Non-redundant uniquely mapped reads","Percentage of reads in chrM","Percentage of reads in chrX","Percentage of reads in chrY","Useful reads","Useful single ends")
	report=append(report,list(Library.size=library))
	} else {
		refer=c(refer,paste(round(mean(ref.chr[20,]/ref.map[,6]),4)*100,'%(SD:',round(sd(ref.chr[20,]/ref.map[,6]),4)*100,'%)',sep=''))
		refer=c(refer,paste(round(mean(ref.map[,8])),'(SD:',round(sd(ref.map[,8])),')',sep=''))
		refer=c(refer,paste(round(mean(ref.useful[,5])),'(SD:',round(sd(ref.useful[,5])),')',sep=''))

		samples=c(samples,paste(round(chr[which(chr[,1]=='chrM'),3]/map[,6],4)*100,'%',sep=''),map[,8],useful[,5])
		library=data.frame(samples,refer)
		colnames(library)=c("Sample","ENCODE PE")
		rownames(library)=c("Total reads","Mapped reads","Non-redundant uniquely mapped reads","Percentage of reads in chrM","Useful reads","Useful single ends")
		report=append(report,list(Library.size=library))
	}

refer=paste(round(mean(ref.dedup[which(ref.dedup$class=='ENCODE PE'),2]),2),'%(SD:',round(sd(ref.dedup[which(ref.dedup$class=='ENCODE PE'),2]),2),'%)',sep='')
refer=c(refer,paste(round(mean(1-ref.map$nodup_ratio),4)*100,'%(SD:',round(sd(1-ref.map$nodup_ratio),4)*100,'%)',sep=''))
samples=c(paste(round(ref.dedup[which(ref.dedup$class=='Sample'),2],2),'%',sep=''),paste((1-map$nodup_ratio)*100,'%',sep=''))
library=data.frame(samples,refer)
colnames(library)=c("Sample","ENCODE PE")
rownames(library)=c("Before alignment library duplicates percentage","After alignment PCR duplicates percentage")
report=append(report,list(Library.complexity=library))

ref.peak=read.table(paste(refpath,'merged_saturation_collection.txt',sep='/'),header=T,sep='\t')
ref.peak=ref.peak[,-1]
ref.peak=ref.peak[11,seq(2,161,3)]
refer=as.numeric(ref.peak[1,])
refer=paste(round(mean(refer)),'(SD:',round(sd(refer)),')',sep='')
refer=c(paste(round(mean(ref.useful[,4]),2)*100,'%(SD:',round(sd(ref.useful[,4]),2)*100,'%)',sep=''),refer)
refer=c(refer,paste(round(mean(ref.map$rup_ratio),2),'%(SD:',round(sd(ref.map$rup_ratio),2),'%)',sep=''))
refer=c(refer,paste(round(mean(ref.enrich2[which(ref.enrich2$class=='ENCODE PE'),1]),2),'(SD:',round(sd(ref.enrich2[which(ref.enrich2$class=='ENCODE PE'),1]),2),')',sep=''))
refer=c(refer,paste(round(mean(ref.enrich3[which(ref.enrich3$class=='ENCODE PE'),1]),2),'(SD:',round(sd(ref.enrich3[which(ref.enrich3$class=='ENCODE PE'),1]),2),')',sep=''))
ref.dicho=read.table(paste(refpath,'merged_bg_dichoto.txt',sep='/'),header=T,sep='\t')
ref.dicho=data.frame(ref.dicho[2:4],class="ENCODE PE")
dicho=read.table(paste("dichoto_bg",name,sep="_"),sep='\t')
dicho=data.frame(dicho,class="Sample")
colnames(dicho)=colnames(ref.dicho)
ref.dicho=rbind(ref.dicho,dicho)
colnames(ref.dicho)[1:3]=c("RPKM smaller than 0.15","RPKM smaller than 0.3","RPKM larger than 0.3")
refer=c(refer,paste(round(mean(ref.dicho[,2]),2),'%(SD:',round(sd(ref.dicho[,2]),2),'%)',sep=''))

samples=c(paste(useful[,4]*100,'%',sep=''),as.numeric(saturate[11,2]))
samples=c(samples,paste(round(map$rup_ratio,2),'%',sep=''),round(ref.enrich2[which(ref.enrich2$class=='Sample'),1],2),round(ref.enrich3[which(ref.enrich3$class=='Sample'),1],2))
samples=c(samples,paste(round(dicho[,2],2),'%',sep=''))

library=data.frame(samples,refer)
colnames(library)=c("Sample","ENCODE PE")
rownames(library)=c("Useful reads ratio","Number of peaks","Reads under peaks ratio","Enrichment ratio in coding promoter regions","Normalized enrichment ratio","Percentage of background RPKM smaller than 0.3777")
report=append(report,list(Enrichment=library))

name=args[6]
capture.output(print(report),file=paste(name,"report.txt",sep='_'))

# json file generation
part1=data.frame("V1",genome,data_type)
colnames(part1)=c("pipe version","genome","read type")
file=list(`data information`=part1)

part2=data.frame("cutadapt","1.12",as.numeric(args[9]),"FastQC","0.11.5")
colnames(part2)=c("program1","program1 version","Removed reads by cutadapt","program2","program2 version")
file=append(file,list(`pre alignment stats`=part2))

part3=data.frame("bwa","0.7.16a","bwa men","methylQA","0.1.9","methylQA atac",map$total,map$mapped,map[,5],map[,6],map[,8],useful[,5])
colnames(part3)=c("alignment program","alignment program version","alignment program parameters","post alignment program","post alignment program version","post alignment program parameters","total reads","mapped reads","uniquely mapped reads","non-redundant mapped reads","useful reads","useful single ends")
file=append(file,list(`mapping stats`=part3))

part9=data.frame(round(ref.dedup[which(ref.dedup$class=='Sample'),2],2)/100,(1-map$nodup_ratio))
colnames(part9)=c("before alignment library duplicates percentage","after alignment PCR duplicates percentage")
file=append(file,list(`library complexity`=part9))

part4=data.frame("plot3.1_insertion_size.png")
colnames(part4)="plot_url"
file=append(file,list(`insert size ditribution`=part4))

autosome=chr[which(!chr$V1%in%c("chrM","chrX","chrY")),c(1,3)]
autosome$V3=round(autosome$V3/map[,6],4)
autosome$V1=paste("@",autosome$V1,"@",sep="")
autosome=paste(autosome$V1,autosome$V3,sep=": ")
autosome=paste("!",paste(autosome,sep="",collapse=", "),"!",sep="")

if(genome!="danRer10") {
	part5=data.frame(round(chr[which(chr[,1]=='chrM'),3]/map[,6],4),round(chr[which(chr[,1]=='chrX'),3]/map[,6],4),round(chr[which(chr[,1]=='chrY'),3]/map[,6],4),autosome)
	colnames(part5)=c("chrM reads percentage","chrX reads percentage","chrY reads percentage","autosome")
	file=append(file,list(`mapping distribution`=part5))
} else {
	part5=data.frame(round(chr[which(chr[,1]=='chrM'),3]/map[,6],4),autosome)
	colnames(part5)=c("chrM reads percentage","autosome")
	file=append(file,list(`mapping distribution`=part5))
}

part6=data.frame("macs2","--keep-dup 1000 --nomodel --shift 0 --extsize 150","qvaule",0.01,map$rup_ratio,map[,11])
colnames(part6)=c("peak calling software","peak calling parameters","peak threshold parameter","peak threshold","reads percentage under peaks","reads number under peaks")
file=append(file,list(`peak analysis`=part6))

part7=data.frame("plot4.5_saturation.png")
colnames(part7)="plot_url"
file=append(file,list(`saturation`=part7))

part8=data.frame(round(ref.enrich2[which(ref.enrich2$class=='Sample'),1],2),round(ref.enrich3[which(ref.enrich3$class=='Sample'),1],2),round(dicho[,2],2)/100)
colnames(part8)=c("enrichment ratio in coding promoter regions","normalized enrichment ratio","percentage of background RPKM smaller than 0.3777")
file=append(file,list(`enrichment`=part8))

part10=data.frame(paste("idr_plot_",name,".ps",sep=""))
colnames(part10)="plot_url"
file=append(file,list(`idr`=part10))

file=list(name=file)
names(file)=name

test=try(library(jsonlite),silent=T)

capture.output(toJSON(file,pretty=T),file=paste(name,"report.json",sep='_'))


