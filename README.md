# Atac-seq_Quality_Control_pipe
This is for the QC matrix construction, data analysis and visualization for Atac-seq data.  
Current version: `V1.1a`  
Current Docker version:   
zhanglab/atac-seq:full `sha256:816fbf3b0e6e3d853b96d5490b0b434f58fb552645d6cfc118e3ff0e2feba8f6`  
zhanglab/atac-seq:mm10 `sha256:521b0d6dad05a0bbf669d07c9d75be4b2d13106c0ba77330593fb9dfb943e42f`  

Advisor: Bo Zhang 
Contributor: Cheng Lyu and Shaopeng Liu  

For any question, please contact Wustl.Zhanglab@gmail.com  


## Usage:  
Singularity 2-step solution (easiest way)  

Step1. download singularity container (you only need download the containcer for once, then you can use them directly):  
#### Please chooice one of them:
```bash
# 1. download from shub (now support mm10 only):  
singularity pull -n zlab_atac.simg shub://Wustl-Zhanglab/ATAC-seq_QC_analysis:mm10  

# 2. download from local server (support mm10, mm9, hg19, hg38 and danRer10):  
wget -O zlab_atac.simg http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/simg/zlab_atac_v1.1a_full.simg  
```

Step2. process data by the singularity image: 
#### Please run at same directory with your data  
```bash
singularity run -B ./:/scratch zlab_atac.simg  -r <SE/PE> -g <mm10/mm9/hg19/hg38/danRer10>  -o <read_file1>  -p <read_file2>  
```

That's it!

#parameters:  
`-r`: SE for single-end, PE for paired-end  
`-g`: genome reference  
`-o`: reads file 1 or the SE reads file, must be ended by .fastq or .fastq.gz or .sra (for both SE and PE)  
`-p`: reads file 2 if this is for PE data, must be ended by .fastq or .fastq.gz  

e.g:
a) mm10 SE data A.fastq  
```bash
singularity run -B ./:/scratch zlab_atac.simg  -r SE -g mm10 -o A.fastq  
```
b) hg38 PE data B_1.fastq B_2.fastq  
```bash
singularity run -B ./:/scratch zlab_atac.simg  -r PE -g hg38 -o B_1.fastq  -p B_2.fastq  
```
c) danRer10 PE data in sra file C.sra  
```bash
singularity run -B ./:/scratch zlab_atac.simg  -r PE -g danRer10 -o C.sra  
```




