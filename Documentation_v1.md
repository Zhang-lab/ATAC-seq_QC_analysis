# Documentation v1
ATAC-seq quality control matrix for Bo Zhang's lab  
Last edit: 10/24/2017  
shaopeng.liu@wustl.edu  											   

## Before use the pipe:  
1, please download [**code_collection**](https://github.com/ShaopengLiu1/Atac-seq_Quality_Control_pipe/tree/master/code_collection).  
2, Please check the [**qc_pipe_source.sh**](https://github.com/ShaopengLiu1/Atac-seq_Quality_Control_pipe/blob/master/code_collection/qc_pipe_source.sh) file to make sure all necessary support files are correctly connected.  
There are 4 species, with an additional whatever-youlike feature. After downloading that file, please revise the link to the resource file in your local server. It's okay to change only 1 of them, as long as you don't need the rest.  
3, to use it, call the "atac_pipe_v1.sh" by bash with proper parameters. Nohup is highly recommended due to the long processing time (approximately 4 hours, the speed limiting process is cutadapt and BWA alignment)

## Pipe Usage:  
user@domain: nohup bash path_to_pipe/pipe.sh  -g  <mm10/mm9/hg38/hg19/personalize>  -r <PE/SE>  -o read_file1  -p read_file2 (if PE file)  &
Optional parameter:   -t <threads>  -m <marker>  -h for help  

## Caveats in method selection
### 1, reads distribution count in each chromosome  
We do **NOT** use samtools index directly, because "bwa" would assign reads of which has mapQ=0 to chr1, the results are not accurate. What we perform is to remove those reads first and count directly by "uniq" command. (This part is consistent with samtools index results)  

### 2, generate random regions from genome  
We do **NOT** use "bedtools random" because it's hard to determine how many regions we need. Especially for Zebra fish data, those peaks(extended to 10kb at each side) would cover most of genome so it's very hard to get enough hits. So we simply shutter the whole genome and take all regions.




## Data Processing:  
### Step1, Pre-alignment   
#### 1.1, Trimming by cutadapt  
Tool: cutadapt v1.12  
input: fastq file  
output: trimmed fastq file  
commands:   
	for PE data: cutadapt -a adapter1 -A adapter2 -q 10 --minimum-length 36  -o $Trimed_output_1  -p $Trimed_output_2  $input_1 $input_2  
	for SE data: cutadapt -a adapter1 -q 10 --minimum-length 36  -o $Trimed_output   $input    
(Default adapter for Atac-seq data: CTGTCTCTTATACACATCT)  
QC to report: cutadapt log  
  
#### 1.2, Fastq file quality assessment  
tool: FastQC v0.11.5  
input: trimmed fastq file  
output: fastqc results  
commands: fastqc -t $threads -o .  $trimed_fastq  
QC to report: fastqc results  


### Step2, reads alignment and reads distribution  
tool:   
  1, bwa v0.7.12  
  2, samtools   
  3, methylQA v0.1.9  
input: trimmed fastq file  
output: aligned bam file  
command:  
bwa mem -t $threads  $mm10_ref_genome.fa  $trimmed.fastq | samtools view -bS - | samtools sort - -O 'bam' -o  $aln.bam -T temp_aln  
QC to report: mapped reads distribution in each chromosome  


### Step3, peak call  
tool:   
  1, methylQA v0.1.9  
  2, macs2 v2.1  
input:   
  1, alined file as $aln.bam  
  2, effect reads file from 1 as $reads.open.bed  
output: peak file  
commands:   
  1, methylQA atac $mm10_ref_genome_size   $aln.bam   
  2, macs2 callpeak -t $reads.open.bed  -g mm -q 0.01 -n peakcall_    --keep-dup 1000 --nomodel --shift 0 --extsize 150  
QC to report:   
	1, mapping status  
	2, insertion size distribution  
note: methylQA atac would do a PE->SE transformation based on insertion site for PE data.  


### Step4, Post alignment analysis:  
#### 4.1, Reads under peak ratio  
tool: bash
input: reads file as $reads.open.bed, and peak file as $peak  
output: reads under peak ratio  
code: rup_number=`intersectBed -a $bed -b $peak -u -f 0.5 | wc -l` 
QC to report: reads under peak ratio  

#### 4.2, Calculate enrichment ration of top 20k peaks against overall background  
tool: bash  
input: reads enrichment in coding promoters, and normalized enrichment for all peak region
output: enrichment ratio file  
QC to report: enrichment ratio for top 20k peak   
  
#### 4.3, PBC calculation (PCR bottlenecking coefficiency)  
tool: bash  
input: mapped reads file   
output: PBC 1  
QC to report: PBC 1  
  
#### 4.4, IDR calculation  
tool: bash and R  
See: [**ENCODE Introduction of IDR**](https://sites.google.com/site/anshulkundaje/projects/idr#TOC-Intuitive-Explanation-of-IDR-and-IDR-plots) for more details  
input: chrom_size file, and pseudo replicates peak file  
ouput:  
  1, IDR peaks   
  2, plot of peak consistency  
 QC to report: IDR results collection  
 
 #### 4.5, Saturation analysis  
 tool: macs2 v2.1  
 input: effect reads  
 commands:  
  1, randomly sampling reads file from 5%, 10%, 20% to 90% of reads file  
  2, call peak and calculate the peaks in sub-sampling file  
 output: saturation result  
 QC to report: saturation plot and table  
  
 #### 4.6, Calculate background  
 tool:   
  1, Python  
  2, intersectBed  
 input: random regions from genome that are at least 10kb from any peak  
 command:  
  1, generate random regions from genome  
  2, remove those with peak in 10kb range  
  3, calculate RPKM for those region  
 output: background record  
 QC to report: background plot  
   
 #### 4.7, plot all results with reference dataset collection (currently ENCODE mm10 PE data)  
 tool: R  
 QC to report: all plots  
 
 
 
