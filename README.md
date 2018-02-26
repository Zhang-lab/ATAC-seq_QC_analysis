# Atac-seq_Quality_Control_pipe
This is for the QC matrix construction and analysis for Atac-seq data.  
Advisor: Bo Zhang  
Contributor: Cheng Lyu and Shaopeng Liu  

## Usage:  
Singularity solution (easiest way)  
1, download singularity container (you only need download the containcer for once, then you can use them directly):  
```bash
singularity pull -n zlab_atac.simg shub://ShaopengLiu1/Zhanglab_ATAC-seq_analysis:mm10  
```

2, process data by the singularity image:  
```bash
singularity exec zlab_atac.simg  -r <SE/PE> -g <mm10/hg38/danRer10>  -o <read_file1>  -p <read_file2>  
```

That's it!

#parameters:  
-r: SE for single-end, PE for paired-end  
-g: genome reference  
-o: reads file 1 or the SE reads file, must be ended by .fastq or .fastq.gz or .sra (for both SE and PE)  
-p: reads file 2 if this is for PE data, must be ended by .fastq or .fastq.gz  

e.g:
a) mm10 SE data A.fastq  
```bash
singularity exec zlab_atac.simg  -r SE -g mm10 -o A.fastq  
```
b) hg38 PE data B_1.fastq B_2.fastq  
```bash
singularity exec zlab_atac.simg  -r PE -g hg38 -o B_1.fastq  -p B_2.fastq  
```
c) danRer10 PE data in sra file C.sra  
```bash
singularity exec zlab_atac.simg  -r PE -g danRer10 -o C.sra  


For any question, please contact Wustl.Zhanglab@gmail.com  



