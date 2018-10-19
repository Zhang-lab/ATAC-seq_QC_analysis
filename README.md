# Atac-seq Integrative Analysis Pipeline  
This is for the QC matrix construction, data analysis and visualization for ATAC-seq data.  
Current version: `IAP_v1.00`   

Advisor: Bo Zhang  
Contributor: Cheng Lyu and Shaopeng Liu  

For any question, please contact Wustl.Zhanglab@gmail.com  


## Usage:  
Singularity 2-step solution (easiest way)  

Step1. download singularity images and reference files (you only need download them **ONCE**, then you can use them directly), if there is any update, you may need to download a new image, but reference files are usually NOT changed:  
####  
1. find and download the image: **[ click here ](http://brc.wustl.edu/SPACE/shaopengliu/Singularity_image/atac-seq/)**, right click to copy the link, and download by wget command. e.g:  
`wget http://brc.wustl.edu/SPACE/shaopengliu/Singularity_image/atac-seq/ATAC_IAP_v1.00.simg`  
2. then **[ click here ](http://brc.wustl.edu/SPACE/shaopengliu/Singularity_image/atac-seq/ref_file/)** to find your interested genome, for now we have mm9/10, hg19/38, danRer10/11, rn6 and dm6. Use the similar way to download them.  

Step2. process data by the singularity image: 
#### Please run the cmd on the same directory of your data, if your data is on /home/example, then you may need `cd /home/example` first. The location of image and reference files is up to you.    
```bash
singularity run -B ./:/process -B <path-to-parent-folder-of-ref-file>:/atac_seq/Resource/Genome  <path-to-downloaded-image> -r <SE/PE> -g <mm10/mm10/hg38 etc.>  -o <read_file1>  -p <read_file2>  
```
It may looks a little confusing at first time, but when you get familier with Singularity they will be friendly :)  
For example, if  
a) you download the image on /home/image/ATAC_IAP_v1.00.simg  
b) the reference file on /home/src/mm10  
c) and your data is read1.fastq.gz and read2.fastq.gz on folder /home/data  

Then you need to:
1. `cd /home/data` 
2. `singularity run -B ./:/process -B /home/src:/atac_seq/Resource/Genome  /home/image/ATAC_IAP_v1.00.simg  -r PE -g mm10 -o read1.fastq.gz -p read2.fastq.gz`  

**explaination**:  
The cmd is in this manner: `singularity run <options> <singularity_image_to_run> <pipeline_parameters>`  

**soft link introduction**:
If you want to use soft link, which is much more friendly when you have a lot of data, of the data. You will only need to add one bind option for singularity, which is `-B <full-path-of-original-position>:<full-path-of-original-position>`  
For example, I want to soft link my data from /scratch to run on my own folder /home/example:  
1. ln -s /scrach/mydata.fastq.gz /home/example; **Please make sure you use the absolute path**  
2. cd /home/example  
3. `singularity run -B ./:/process -B /home/src:/atac_seq/Resource/Genome -B /scratch:/scratch  /home/image/ATAC_IAP_v1.00.simg  -r PE -g mm10 -o read1.fastq.gz -p read2.fastq.gz`  

#parameters:  
`-h`: help information  
`-r`: SE for single-end, PE for paired-end  
`-g`: genome reference, one simg is designed for ONLY one species due to the file size. For now the supported genoms are: <mm10/mm9/hg19/hg38/danRer10/danRer11/rn6/dm6>  
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





