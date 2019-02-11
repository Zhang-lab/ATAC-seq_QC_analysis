# Documentation v1.00
ATAC-seq quality control and analysis pipeline for Bo Zhang's lab  
Last edit: 02/11/2019  
For any question please contact: shaopeng.liu@wustl.edu  

**Outline**  
I, Output example and annotation  
II, Terms  
III, Data processing details  

## I, Output example  
After running the pipeline, there will be a folder called **Processed_${name}**, all intermediate files and final output files are stored there. And it looks like this (using an ENCODE data for example):  
![output image](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/AIAP_output_example.png)  

There should be 15 files in total, and please feel free to **[ click here ](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/Processed_mm10_liver_embryo_d11.5_PE_bio1_ENCLB441LCB_1/)** to explore them.  

Regarding each file (${name} refers to the file name, at here it's "mm10_liver_embryo_d11.5_PE_bio1_ENCLB441LCB_1"):  
File name | Content  
--------- | -------  
*QC_ATAC_data_collection_${name}* | stores all intermediate output files from each step, and they are merged into one output json file and txt file in the pipeline  
*QC_${name}.json* | record all QC information of current data and ENCODE reference data (by pipe v3.1, will update soon)  




   
2, QC_pipe_processing.log: store the status of each step, and warning messages if any  
3, QC_data_collection_${file}.result: this table is inside the QC data collection folder, it's a one line table with same information as the json file but easier to collect in batch for review on server  
4, Single output: for each step, the intermediate files are kept  

## II, Term and definition  
1, coding promoter region definition:  
2kb region (1kb up and 1kb down) of the transcription start site, this is a rough estimate and doesn't account multiple starting point situation.  

2, enrichment in coding promoter region:  
coding enrichment = ( reads in promoter / promoter length)  /  (total reads / genome size)  

3, insertion length (from output of methylQA atac):  
The distance between 2 insertion site from the ATAC-seq reads.  

4, peak length (from output of macs2):  
The distance between start point and end point of a peak.  

5, unique chrM ratio:  
The percentage of uniquely mapped chrM reads in all uniquely mapped reads. It's used to measure the sample quality.  

6, mapping status:  
> total reads: the raw reads of the input fastq file  
> written reads: the output reads from trimming  
> mapped reads: the reads that can be mapped to reference genome  
> uniquely mapped reads: the reads that can be mapped uniquely (mapQ > 10)  
> useful reads: non-redundant uniquely mapped reads  
> useful single ends: in ATAC-seq, what we care about is the insertion rather than the fregment itself. So we do a transformation on the reads (please see methylQA atac for more information) to focusing on the insertion point only. The reads after this transformation is called "useful single ends".  

7, reads under peak:  
Those reads whose center point located inside a peak region.  

8, subsample 10 million enrichment:  
This is how we estimate the enrichment of reads for each single data.  
We will subsample the useful single ends down to 10M, and calculate based on the subset of data.  
sub 10M enrichment =  
(($rupn+10000000*$peak_length / $genome_size) / $peak_length)   /    (20M / ($genome_size-$peak_length))  

9, background RPKM:  
> Random sample 500bp regions from genome
> Keep those region that are far from any know peaks (distance > 10kb)  
> Calculate the RPKM for those kept regions  


## III, Data Processing  
### Caveats in method selection
1, reads distribution count in each chromosome  
We do **NOT** use samtools index directly, because "BWA" would assign reads of which has mapQ=0 to chr1, the results are not accurate. What we perform is to remove those reads first and count directly by "uniq" command. (This part is consistent with samtools index results)   
2, generate random regions from genome to calculate background  
We do **NOT** use "bedtools random" because it's hard to determine how many regions we need. Especially for Zebra fish data, those peaks(extended to 10kb at each side) would cover most of genome so it's very hard to get enough hits. So we simply shutter the whole genome and take all regions.  

### Step1, Pre-alignment   
#### 1.1, Trimming by cutadapt  
Tool: cutadapt v1.16  
input: fastq file  
output: trimmed fastq file  
commands:   
for PE data: 
> cutadapt -a adapter1 -A adapter2 -q 10 --minimum-length 36  -o $Trimed_output_1  -p $Trimed_output_2  $input_1 $input_2  
for SE data: 
> cutadapt -a adapter1 -q 10 --minimum-length 36  -o $Trimed_output   $input    
(Default adapter for Atac-seq data: CTGTCTCTTATACACATCT)  
QC to report: cutadapt log  
  
#### 1.2, Fastq file quality assessment  
tool: FastQC v0.11.7  
input: trimmed fastq file  
output: fastqc results  
commands: 
> fastqc -t $threads -o .  $trimed_fastq  
QC to report: fastqc results  


### Step2, reads alignment and reads distribution  
#### 2.1, BWA MEM aligner  
tool: bwa v0.7.12  
input: trimed fastq file  
commands:  
> bwa mem $bwa_ref 'step1.1_trimed_'$name*'.fastq' | $samtools view -bS - | $samtools sort - -O 'bam' -o 'step2.1_trimed_'$name'.bam' -T temp_aln  
output: aligned reads in bam file  

#### 2.2, reads distribution  
tool: samtools 1.3.1, methylQA v0.2.1    
input: aligned bam file  
output: reads distribution count  
QC to report: mapped reads distribution in each chromosome  

#### 2.3, library complexity estimate  
tool: preseq 2.0.0 (from subread)  
input: aligned bam file  
commands:  
> preseq lc_extrap -o 'step2.3_yield_'$name'.result' -B  'step2.1_trimed_'$name'.bam'  
output: library complexity estimate  
QC to report: library complexity estimate  

### Step3, peak call  
tool:   
  1, methylQA v0.1.9  
  2, macs2 v2.1  
input:   
  1, alined file as $aln.bam  
  2, effect reads file from 1 as $reads.open.bed  
output: peak file  
commands:   
> methylQA atac $mm10_ref_genome_size   $aln.bam   
> macs2 callpeak -t $reads.open.bed  -g mm -q 0.01 -n peakcall_    --keep-dup 1000 --nomodel --shift 0 --extsize 150  
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
  
#### 4.3 (deleted), PBC calculation (PCR bottlenecking coefficiency)  
tool: bash  
input: mapped reads file   
output: PBC 1  
QC to report: PBC 1  
  
#### 4.4, Saturation analysis  
tool: macs2 v2.1  
input: effect reads  
commands:  
> randomly sampling reads file from 5%, 10%, 20% to 90% of reads file  
> call peak and calculate the peaks in sub-sampling file  
output: saturation result  
QC to report: saturation plot and table  
  
#### 4.5, Calculate background  
tool:   
  1, Python  
  2, intersectBed  
input: random regions from genome that are at least 10kb from any peak  
command:  
> generate random regions from genome  
> remove those with peak in 10kb range  
> calculate RPKM for those region  
output: background record  
QC to report: background plot  
   
#### 4.6, plot all results with reference dataset collection (currently ENCODE mm10 PE data)  
tool: R  
QC to report: all plots and report file 
 
 
 
