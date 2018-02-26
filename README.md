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



For any question, please contact Wustl.Zhanglab@gmail.com  



