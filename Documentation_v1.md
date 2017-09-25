# Documentation v1
ATAC-seq quality control matrix for Bo Zhang's lab  
Last edit: 09/25/2017  
shaopeng.liu@wustl.edu  											   

## Before use the pipe:  
Please check the **qc_pipe_source.sh** file to make sure all necessary support files are correctly connected.  

## Pipe Usage:  
user@domain: path_to_pipe/pipe.sh  -g <genome>  -r <PE/SE>  -o read_file1  -p read_file2 (if PE file)  
Optional parameter:   -t <threads>  -m <marker>  -h for help  

## Data Processing:
### Step1, Pre-alignment Trimming.
Tool: cutadapt v1.12
