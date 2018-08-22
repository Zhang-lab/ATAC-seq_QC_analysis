# Atac-seq_Quality_Control_pipe
This is for the QC matrix construction, data analysis and visualization for ATAC-seq data.  
Current version: `v4`   

Advisor: Bo Zhang  
Contributor: Cheng Lyu and Shaopeng Liu  

For any question, please contact Wustl.Zhanglab@gmail.com  


## Usage:  
Singularity 2-step solution (easiest way)  

Step1. download singularity container (you only need download the containcer for once, then you can use them directly):  
####  
```bash
# download image from local server:  
wget http://brc.wustl.edu/SPACE/shaopengliu/Singularity_image/atac-seq/ATAC_seq_mm10_v4.simg  
```

Step2. process data by the singularity image: 
#### Please run at same directory with your data OR the soft link of your data    
```bash
singularity run -H ./:/scratch <path-to-image> -r <SE/PE> -g <mm10>  -o <read_file1>  -p <read_file2>  
```

That's it!

#parameters:  
`-h`: help information  
`-r`: SE for single-end, PE for paired-end  
`-g`: genome reference, one simg is designed for ONLY one species due to the file size. For now the supported genoms are: <mm10/mm9/hg19/hg38/danRer10> (only mm10 in the example).  
`-o`: reads file 1 or the SE reads file, must be ended by .fastq or .fastq.gz or .sra (for both SE and PE)  
`-p`: reads file 2 if input PE data, must be ended by .fastq or .fastq.gz  
`-c`: (optional) specify read length minimum cutoff for methylQA filtering, default 38  
`-t`: (optional) specify number of threads to use, default 24  
`-i`: (optional) insertion free region finding parameters used by Wellington Algorithm (Jason Piper etc. 2013), see documentation for more details.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If you don NOT want to run IFR finding step, please just ignore the -i option; however IFR finding will use default parameters only if -i specified as 0:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;min_lfp=5  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;max_lfp=15  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;step_lfp=2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;min_lsh=50  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;max_lsh=200  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;step_lsh=20  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;method=BH  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;p_cutoff=0.05  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If you want to specify your own parameter, please make sure they are in the same order and seperated by comma  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example: -i 5,15,2,50,200,20,BH,0.05  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You can check the pipe log file for the parameters used by IFR code  

e.g:
a) mm10 SE data A.fastq  
```bash
singularity run -H ./:/scratch ATAC-seq_v3.1_mm10.simg  -r SE -g mm10 -o A.fastq  
```
b) hg38 PE data B_1.fastq B_2.fastq  
```bash
singularity run -H ./:/scratch ATAC-seq_v3.1_hg38.simg  -r PE -g hg38 -o B_1.fastq  -p B_2.fastq  
```
c) danRer10 PE data in sra file C.sra  
```bash
singularity run -H ./:/scratch ATAC-seq_v3.1_danRer10.simg  -r PE -g danRer10 -o C.sra  
```




