#!/usr/bin/env Rscript
args=commandArgs()

# R script using Wellington algorithm to calculate the p-value of insertion free region #
#########################################################################################
library(data.table)
library(Rcpp)

# cpp implement of `outer(vec1,vec2,"<=")` in R
cppFunction("
  NumericMatrix cNoLarger(NumericVector vec1, NumericVector vec2) {
    int n1 = vec1.size();
    int n2 = vec2.size();
    NumericMatrix out(n1,n2);
    for(int i = 0; i < n1; i++) {
      for(int j =0 ; j < n2; j++) {
        if(vec1[i] <= vec2[j]) {
          out(i,j) = 1;
        }
      }
    }
    return out;
  }
")

# cpp implement of `outer(vec1,vec2,">=")` in R
cppFunction("
  NumericMatrix cNoSmaller(NumericVector vec1, NumericVector vec2) {
    int n1 = vec1.size();
    int n2 = vec2.size();
    NumericMatrix out(n1,n2);
    for(int i = 0; i < n1; i++) {
      for(int j =0 ; j < n2; j++) {
        if(vec1[i] >= vec2[j]) {
          out(i,j) = 1;
        }
      }
    }
    return out;
  }
")

# cpp implement of `X%*%y` in R
cppFunction("
  NumericVector cMatVecMult(NumericMatrix X, NumericVector y){
    int ncol = X.ncol();
    int nrow = X.nrow();
    NumericVector out(nrow);
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        out[i] += X(i,j) * y[j];
      }
    }
    return out;
  }
")

# find one insertion-free given a specific region range
find_motif=function(start_peak,end_peak,bed,min_lfp,max_lfp,step_lfp,min_lsh,max_lsh,step_lsh,user_method,p_cutoff) {
  # create all possible combination of center `c`, half footprint length `l` and shoulder length `s`
  comb=expand.grid(seq(start_peak+min_lsh+max_lfp,end_peak-min_lsh-max_lfp,1),seq(min_lfp,max_lfp,step_lfp),seq(min_lsh,max_lsh,step_lsh))
  minus=comb[,1]-comb[,2] # c-l
  plus=comb[,1]+comb[,2] # c+l
  max=pmax(comb[,1]-comb[,2]-comb[,3],start_peak) # max(c-l-s,start_peak)
  min=pmin(comb[,1]+comb[,2]+comb[,3],end_peak) # min(c+l+s,end_peak)
  p=comb[,2]/(comb[,2]+comb[,3])
  
  # create boolean matrix for better speed using matrix multiplication
  bool1=cNoLarger(minus,bed$V2) # (c-l)<=bed$V2
  bool2=cNoSmaller(plus,bed$V2) # (c+l)>=bed$V2
  bool3=cNoLarger(max,bed$V2) # max(c-l-s,start_peak)<=bed$V2
  bool4=cNoSmaller(min,bed$V2) # min(c+l+s,end_peak)>=bed$V2
  
  # count fp,sh_d,sh_u for each combination of (c,l,s)
  fp=cMatVecMult((bool1*bool2),bed$V4)
  sum=fp+cMatVecMult((bool3*(1-bool1)),bed$V4)+cMatVecMult(((1-bool2)*bool4),bed$V4)
  binomial=as.data.frame(cbind(fp,sum,p))
  pval=with(binomial,pbinom(fp,sum,p))
  
  # find the corresponding (l,s) of each `c` with the smallest p-value
  re=cbind(comb,pval)
  re=data.table(re)
  re=re[,.SD[which.min(pval)],by=Var1]
  
  # adjust p-value for given method
  re$padj=p.adjust(re$pval,method=user_method)
  smallest=re[order(re$pval)[1],]
  
  if(smallest$padj<p_cutoff) {
    return(smallest[1,c(1,2,5)])
  } else {
    return(0)
  }
}

find_IF=function(count,min_lfp=5,max_lfp=15,step_lfp=1,min_lsh=50,max_lsh=200,step_lsh=2,user_method="BH",p_cutoff=0.05) {
  IF_regions=c()
  peak_index=names(table(count[,8]))
  # scan all peaks
  for(peak in peak_index) {
    bed=count[count[,8]==peak,]
    region_list=data.frame(start=bed[1,6],end=bed[1,7])
    bed0=bed[,1:4]
    chr=bed[1,1]
    # find all possible motif region wihtin a peak
    while(dim(region_list)[1]!=0) {
      start_peak=region_list[1,1]
      end_peak=region_list[1,2]
      bed=bed0[bed0$V2>=start_peak & bed0$V3<=end_peak,]
      if(start_peak+min_lsh+max_lfp<=end_peak-min_lsh-max_lfp) {
        motif=as.data.frame(find_motif(start_peak,end_peak,bed,min_lfp,max_lfp,step_lfp,min_lsh,max_lsh,step_lsh,user_method,p_cutoff))
        if(length(motif)>1) { 
          IF_regions=rbind(IF_regions,data.frame(chr=chr,start=(motif[1,1]-motif[1,2]),end=(motif[1,1]+motif[1,2]),pvalue=motif[1,3],peak_index=peak))
          # remove found motif from peak region
          region_list=rbind(region_list,data.frame(start=c(region_list[1,1],motif[1,1]+motif[1,2]),end=c(motif[1,1]-motif[1,2],region_list[1,2])))
        }
      }
      region_list=region_list[-1,]
    }
  }
  return(IF_regions)
}

count=read.table(args[6],colClasses=c("character",rep("numeric",3),"character",rep("numeric",2),"character"))

a=find_IF(count,min_lfp=as.numeric(args[8]),max_lfp=as.numeric(args[9]),step_lfp=as.numeric(args[10]),min_lsh=as.numeric(args[11]),max_lsh=as.numeric(args[12]),step_lsh=as.numeric(args[13]),user_method=args[14],p_cutoff=as.numeric(args[15]))

write.table(a,args[7],sep="\t",col.names=F,row.names=F)

