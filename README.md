# Atac-seq_Quality_Control_pipe
This is for the QC matrix construction, data analysis and visualization for Atac-seq data.  
Current version: `V3.1`   

Advisor: Bo Zhang  
Contributor: Cheng Lyu and Shaopeng Liu  

For any question, please contact Wustl.Zhanglab@gmail.com  


## Usage:  
Singularity 2-step solution (easiest way)  

Step1. download singularity container (you only need download the containcer for once, then you can use them directly):  
#### Please chooice one of them:
```bash
# 1. download from local server (mm10 image):  
wget -O zlab_atac.simg http://brc.wustl.edu/SPACE/shaopengliu/Singularity_image/atac_mm10_v3.1.simg  
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
`-t`: (optional) specify threads, default 24  

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




