#!/usr/bin/env Rscript
args=commandArgs()

name=args[6]
name=paste(name,"result",sep=".")

library(ggplot2)
library(cowplot)

# chromosome distribution
chr=read.table(paste("chrom_count",name,sep="_"))
r=as.numeric(rownames(chr[which(chr$V1=='chrM'),]))
chr=chr[c(1:r-1,range(r+1,nrow(chr)),r),]
rownames(chr)=NULL
chr$V1=factor(c(nrow(chr):1),labels=rev(as.character(chr$V1)))
chr$V2=chr$V2/sum(chr$V2)
chr$V3=chr$V3/sum(chr$V3)
colnames(chr)=c("chromosome","Total reads","Non-redundant uniquely mapped reads")
chr=data.frame(chr[1],stack(chr[2:3]))

png("stacked_barplot_chromosome_count.png",height=640,width=6600,res=300)
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
saturate=read.table(paste("saturation",name,sep="_"),header=T)
saturate=saturate[,-2]
colnames(saturate)=c("depth","peak","percentage","marker")

png("line_chart_saturation_peak_percentage.png",height=2300,width=2800,res=300)
ggplot(saturate,aes(x=depth,y=100*percentage))+
	geom_point()+geom_line(size=1)+expand_limits(x=0,y=0)+
	geom_text(aes(label=peak),hjust=0,vjust=1.4,family="Tahoma",fontface=2,size=3)+
	ggtitle("Line chart of recaptured peaks percentage by subsampling original library")+
	scale_y_continuous(name="Percentage of recaptured peaks in original peaks",limits=c(0,100))+
	scale_x_continuous(name="Percentage of original library size",breaks=seq(0,100,by=10))+
	geom_vline(xintercept=50,linetype="dotted",size=1)+
	theme_bw()+theme_classic()+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"))
dev.off()

# peak legth distribution
peakl=read.table(paste("peak_length_distri",name,sep="_"))
peakl=rbind(peakl[which(peakl$V1<1500),],c(1500,sum(peakl[which(peakl$V1>=1500),2])))
peakl$V2=peakl$V2/sum(peakl$V2)

png("density_plot_peak_length.png",height=2000,width=3000,res=300)
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

png("density_plot_insertion_size.png",height=1800,width=2400,res=300)
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
ref.dedup=read.table('/home/chengl/Update/visual_reference/merged_dedup_percentage.txt',header=T) ####
ref.dedup[seq(1,nrow(ref.dedup),2),4]="Fastq1"
ref.dedup[seq(2,nrow(ref.dedup),2),4]="Fastq2"
colnames(ref.dedup)[4]='Reference: ENCODE PE'
dedup=read.table(paste("dedup_percentage",name,sep="_"),header=T)

png("density_plot_dedup.png",height=2000,width=3300,res=300)
ggplot(ref.dedup,aes(x=deduplication_percentage,y=..scaled..,fill=`Reference: ENCODE PE`))+
	geom_density(size=1,adjust=1,alpha=0.5)+
	ggtitle("Density plot of deduplication percentage (Adjust=1)")+
	theme_bw()+theme_classic()+
	scale_y_continuous(name="Density")+
	scale_x_continuous(name="Deduplication percentage",limits=c(0,100))+
	geom_vline(xintercept=dedup[1,2],linetype="dotted",size=1,color="#E69F00")+
	geom_vline(xintercept=dedup[2,2],linetype="dotted",size=1,color="#56B4E9")+
	geom_text(aes(x=dedup[1,2],y=1.2,label=paste("Sample fastq1",round(dedup[1,2],2),sep="\n")),hjust=1.1,size=3.5,family="Tahoma",fontface="bold",color="#E69F00")+
	geom_text(aes(x=dedup[2,2],y=1.1,label=paste("Sample fastq2",round(dedup[2,2],2),sep="\n")),hjust=1.1,size=3.5,family="Tahoma",fontface="bold",color="#56B4E9")+
	scale_fill_manual(values=c("#E69F00","#56B4E9"))+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"),
		legend.text=element_text(size=10,face="bold"),
		legend.title=element_text(size=10,face="bold"),
		legend.position="bottom")
dev.off()

# enrichment
ref.enrich=read.table('/home/chengl/Update/visual_reference/merged_enrichment_ratio.txt',header=T) ####
ref.enrich=ref.enrich[c(4,39),-1]
enrich=read.table(paste("enrichment",name,sep="_"),header=T)
enrich=enrich[,-1]
ref.enrich=data.frame(x=seq(1000,20000,by=1000),upper=as.numeric(ref.enrich[2,]),lower=as.numeric(ref.enrich[1,]),sample=as.numeric(enrich[1,]))

png("line_chart_enrichment_peaks.png",height=2300,width=2600,res=300)
plot=ggplot(ref.enrich,aes(x=x,y=sample))+
	geom_ribbon(aes(ymin=lower,ymax=upper,x=x,fill="Range of ENCODE PE enrichment ratio"),alpha=0.5)+
	geom_point()+geom_line(size=1)+expand_limits(x=0,y=0)+
	ggtitle("Line chart of peaks enrichment")+
	scale_x_continuous(name="Number of Top peaks",limits=c(0,20000))+
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
if(ref.enrich[1,4]<200) {
	plot+scale_y_continuous(name="Enrichment ratio",limits=c(0,200))
} else if(ref.enrich[1,4]<400) {
	plot+scale_y_continuous(name="Enrichment ratio",limits=c(0,400))
} else {
	plot+scale_y_continuous(name="Enrichment ratio")
}
dev.off()

# mapping status
ref.map=read.table('/home/chengl/Update/visual_reference/merged_mapping_status.txt',header=T) ####
ref.map=ref.map[,c(2,7,10,12)]
map=read.table(paste("mapping_status",name,sep="_"),header=T)
map=map[,c(2,7,10,12)]

png("density_plot_noduplication_ratio.png",height=1600,width=2300,res=300)
plot=ggplot(ref.map,aes(x=nodup_ratio,y=..scaled..,fill="Distribution of ENCODE PE noduplication ratio"))+
	geom_density(size=1,adjust=1)+
	ggtitle("Density plot of noduplication ratio (Adjust=1)")+
	theme_bw()+theme_classic()+
	scale_y_continuous(name="Density")+
	geom_vline(xintercept=map[1,3],linetype="dotted",size=1)+
	geom_text(aes(x=map[1,3],y=1.2,label=paste("Sample",map[1,3],sep="\n")),hjust=1.1,size=3.5,family="Tahoma",fontface="bold")+
	scale_fill_manual("Reference density",values="grey")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"),
		legend.text=element_text(size=8,face="bold"),
		legend.title=element_text(size=12,face="bold"),
		legend.position="bottom")
if(map[1,3]>=0.5) {
	plot+scale_x_continuous(name="Noduplication ratio",limits=c(0.5,1))
} else {
	plot+scale_x_continuous(name="Noduplication ratio",limits=c(0,1))
}
dev.off()

png("density_plot_rup_ratio.png",height=1600,width=2300,res=300)
plot=ggplot(ref.map,aes(x=rup_ratio,y=..scaled..,fill="Distribution of ENCODE PE reads under peaks ratio"))+
	geom_density(size=1,adjust=1)+
	ggtitle("Density plot of reads under peaks ratio (Adjust=1)")+
	theme_bw()+theme_classic()+
	scale_y_continuous(name="Density")+
	geom_vline(xintercept=map[1,4],linetype="dotted",size=1)+
	geom_text(aes(x=map[1,4],y=1.2,label=paste("Sample",map[1,4],sep="\n")),hjust=-0.1,size=3.5,family="Tahoma",fontface="bold")+
	scale_fill_manual("Reference density",values="grey")+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"),
		legend.text=element_text(size=8,face="bold"),
		legend.title=element_text(size=12,face="bold"),
		legend.position="bottom")
if(map[1,4]<=60) {
	plot+scale_x_continuous(name="Reads under peaks ratio",limits=c(0,60))
} else {
	plot+scale_x_continuous(name="Reads under peaks ratio",limits=c(0,100))
}
dev.off()

colnames(ref.map)[1:2]=c("Total reads","Non-redundant uniquely mapped reads")
ref.map=stack(ref.map[,1:2])
colnames(ref.map)[2]="Reference: ENCODE PE"

png("density_plot_reads_distribution.png",height=2000,width=3000,res=300)
plot=ggplot(ref.map,aes(x=values,y=..scaled..,fill=`Reference: ENCODE PE`))+
	geom_density(size=1,adjust=1,alpha=0.5)+
	ggtitle("Density plot of reads distribution (Adjust=1)")+
	theme_bw()+theme_classic()+
	scale_y_continuous(name="Density")+
	geom_vline(xintercept=map[1,2],linetype="dotted",size=1,color="#E69F00")+
	geom_vline(xintercept=map[1,1],linetype="dotted",size=1,color="#56B4E9")+
	geom_text(aes(x=map[1,2],y=1.2,label=paste("Sample non-redundant\nuniquely mapped reads",map[1,2],sep="\n")),hjust=1.1,size=3.5,family="Tahoma",fontface="bold",color="#E69F00")+
	geom_text(aes(x=map[1,1],y=1.05,label=paste("Sample total reads",map[1,1],sep="\n")),hjust=1.1,size=3.5,family="Tahoma",fontface="bold",color="#56B4E9")+
	scale_fill_manual(values=c("#E69F00","#56B4E9"))+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"),
		legend.text=element_text(size=10,face="bold"),
		legend.title=element_text(size=10,face="bold"),
		legend.position="bottom")
if(map[1,2]<9e+7) {
	plot+scale_x_continuous(name="Number of reads",limits=c(0,10e+7))
} else {
	plot+scale_x_continuous(name="Number of reads")
}
dev.off()

# promoter
promoter=read.table(paste("promoter_percentage",name,sep="_"),header=T)
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

png("pie_chart_promoter_region.png",height=1500,width=5000,res=300)
options(warn=-1)
plot_grid(p2,p1,ncol=2,nrow=1,rel_widths=c(1,0.88))
options(warn=0)
dev.off()

top=read.table(paste("bin",name,sep="_"))
top=data.frame(rank=as.numeric(top$V1),index=as.numeric(top$V2))

png("promoter_distribution.png",height=1800,width=3200,res=300)
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

ref.background=read.table('/home/chengl/Update/visual_reference/merged_background.txt') ####
ref.background=data.frame(RPKM=ref.background$V6,group=0)
background=read.table(paste("background",name,sep="_"))
background=data.frame(RPKM=background$V6,group=1)
ref.background=rbind(ref.background,background)
ref.background$group=factor(ref.background$group,label=c("ENCODE PE","Sample"))

png("density_plot_background.png",height=2000,width=3300,res=300)
options(warn=-1)
ggplot(ref.background,aes(x=RPKM,fill=group))+
	geom_density(size=1,adjust=1,alpha=0.5)+
	ggtitle("Density plot of background RPKM (Adjust=1)")+
	theme_bw()+theme_classic()+
	scale_y_continuous(name="Density",limits=c(0,10))+
	scale_x_continuous(name="RPKM",limits=c(0,1))+
	geom_vline(xintercept=0.336,linetype="dotted",size=1)+
	geom_text(aes(x=0.336,y=9,label="Theoretical RPKM of background"),hjust=-0.1,size=3.5,family="Tahoma",fontface="bold")+
	scale_fill_manual(values=c("#E69F00","#56B4E9"))+
	theme(plot.title=element_text(size=14,family="Tahoma",face="bold",hjust=0.5),
		text=element_text(size=12,family="Tahoma"),
		axis.title=element_text(face="bold"),
		axis.text.x=element_text(size=10,face="bold"),
		axis.text.y=element_text(size=10,face="bold"),
		legend.text=element_text(size=10,face="bold"),
		legend.title=element_text(size=10,face="bold"),
		legend.position="bottom")
options(warn=0)
dev.off()




