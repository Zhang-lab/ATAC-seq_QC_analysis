# Documentation v1.00
ATAC-seq quality control and analysis pipeline for Bo Zhang's lab  
Last edit: 02/11/2019  
For any question please contact: shaopeng.liu@wustl.edu  

**Outline**  
I, Output example and annotation  
II, Visulization through qATACViewer
III, Terms  
IV, Data processing details  

## I, Output example  
After running the pipeline, there will be a folder called **Processed_${name}**, all intermediate files and final output files are stored there. And it looks like this (using an ENCODE data for example):  

There should be 15 files in total, and please feel free to **[ click here ](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/Processed_mm10_liver_embryo_d11.5_PE_bio1_ENCLB441LCB_1/)** to explore them.  
 
File name | Content
--------- | -------
*QC_ATAC_data_collection_${name}* | stores all intermediate output files from each step, and they are merged into one output json file and txt file in the pipeline
*QC_${name}.json* | record all QC information of current data and ENCODE reference data (by pipe v3.1, will update soon)
*QC_pipe_processing.log* | record the failure status of each step, this signal is captured by $?
*step1.1_${name}_cutadapt_PE.trimlog* | cutadapt trimming report
*step2.1_trimed_${name}.bam* | aligned bam file from **BWA MEM**, keep this one only for backup purpose
*step3.1_normalized_per_10M_${name}.bigWig* | normalized signal per 10M input single ends for visualization purpose
*step3.2_rmbl_${name}.bigWig* | full signal for visualization purpose
*step3.3_rmbl_${name}.open.bed* | bed file after quality filtering on aligned bam file and signal shifting, this is the input for **macs2**
*step3.4_peakcall_${name}_peaks.narrowPeak* | macs2 output, bed file that records peaks called by the software
*step3.4_peakcall_${name}_peaks.xls* | macs2 output, record the command and corresponding output
*step3.4_peakcall_${name}_summits.bed* | macs2 output, bed file that record summits only for each peak
*step4.2_insertion_site_${name}.bedGraph* | bedGraph file that stores all insertion sites
*step4.2_insertion_site_${name}.bigWig* | same to previous file but for visulization use
*step4.6_multiqc_data* | output folder from **multiqc**, please see its help page for more details
*step4.6_multiqc_report.html* | together with previous folder, this html visulize many quality control information

If you happens to open the first folder which record intermediate files, you will see:
![QC image](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/QC_col_folder.png)  

There will be 33 files, but don't worry, all of them are summaized together in the json file. You might be interested some figures stored in the first folder "plots_collection_${name}". If you open it, you will get this.
![plot collection](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/plot_collection.png)  

File name | Content
--------- | -------
*plot2.2_reads_distri_in_chrom.png* | distribution of each chromosome
*plot2.3_yield_distinction.png* | expected library complexity from **preseq**
*plot3.1_insertion_size.png* | smoothed density plot of insertion size (length between 2 insertion points)
*plot3.1_library_reads_distri.png* | library size compared with ENCODE reference
*plot3.3_peak_length.png* | peak length distribution
*plot4.1_RUP.png* | reads under peak ratio compared with ENCODE reference
*plot4.2.2_peaks_enrichment_ratio.png* | enrichment ratio compared with ENCODE reference
*plot4.3_PCR_duplicates_percentage.png* | PCR duplicate ratio compared with ENCODE reference
*plot4.5_saturation.png* | saturation plot
*plot4.6_promoter-peak_count.png* | pie plot of reads and peaks distribution in / out promoter regions
*plot4.6_promoter_distribution_among_peaks.png* | percentage of peaks that cover a promoter in ranked peaks (by qvalue)



## II, Visulization through qATACViewer
To make better visulization, we have also prepared a tool named **qATACViewer**, please **[click here](https://github.com/lidaof/qATACviewer/tree/localjson)** to find its github page. 

**Usage**
1. copy the git repository to local: `git clone -b localjson https://github.com/lidaof/qATACviewer.git` (please makesure to use the `localjson` branch)  
2. update the path in `./qATACViewer/frontend/src/data.json` **(we are improving this part)**    
3. follow the instruction on the the **Readme** of it  



## III, Term and definition  
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


## IV, Data Processing  
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
 
 
 
