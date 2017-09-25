# Documentation v1
ATAC-seq quality control matrix for Bo Zhang's lab  
Last edit: 09/25/2017  
shaopeng.liu@wustl.edu  											   

## Before use the pipe:  
Please check the [**qc_pipe_source.sh**](https://github.com/ShaopengLiu1/Atac-seq_Quality_Control_pipe/blob/master/code_collection/qc_pipe_source.sh) file to make sure all necessary support files are correctly connected.  
There are 4 species, with an additional whatever-youlike feature. After downloading that file, please revise the link to the resource file in your local server. It's okay to change only 1 of them, as long as you don't need the rest.

## Pipe Usage:  
user@domain: path_to_pipe/pipe.sh  -g <genome>  -r <PE/SE>  -o read_file1  -p read_file2 (if PE file)  
Optional parameter:   -t <threads>  -m <marker>  -h for help  

## Data Processing:
### Step1, Pre-alignment Trimming.
Tool: cutadapt v1.12
