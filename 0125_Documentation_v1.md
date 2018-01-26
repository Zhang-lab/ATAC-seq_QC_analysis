# Documentation v1
ATAC-seq quality control matrix for Bo Zhang's lab  
Last edit: 01/25/2018  
For any question please contact: shaopeng.liu@wustl.edu  											   

**Outline**:  
I, Term and definition  
II, Output results  
III, Data processing steps   


## I, Term and definition  
1, coding promoter region definition:  
> Abstract all coding genes from GTF file  
> Merge all transcript region to get a rough estimate on TSS  
> choose 1kb up/downstream of the TSS as the coding promoter region  

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
> mapped reads: the reads that can be mapped to reference genome  
> uniquely mapped reads: the reads that can be mapped uniquely (mapQ > 10)  
> useful reads: non-redundant uniquely mapped reads  
> useful single ends: in ATAC-seq, what we care about is the insertion rather than the fregment itself. So we do a transformation on the reads (please see methylQA atac for more information) to focusing on the insertion point only. The reads after this transformation is called "useful single ends".  

7, subsample 10 million enrichment:  
This is how we estimate the enrichment of reads for each single data.  
We will subsample the useful single ends down to 10M, and calculate based on the subset of data.
sub 10M enrichment =  
(($rupn+10000000*$peak_length / $genome_size) / $peak_length)   /    (20M / ($genome_size-$peak_length))  

8, background RPKM:  
> Random sample 500bp regions from genome
> Keep those region that are far from any know peaks (distance > 10kb)  
> Calculate the RPKM for those kept regions  


## II, Output results  
1, There would be 2 files ended with "report.txt" and "json" that record all related QC information inside.  
2, For every single results in detail, please go to the folder "result_collection".  


## III, Data Processing  
### Caveats in method selection
1, reads distribution count in each chromosome  
We do **NOT** use samtools index directly, because "BWA" would assign reads of which has mapQ=0 to chr1, the results are not accurate. What we perform is to remove those reads first and count directly by "uniq" command. (This part is consistent with samtools index results)   
2, generate random regions from genome to calculate background  
We do **NOT** use "bedtools random" because it's hard to determine how many regions we need. Especially for Zebra fish data, those peaks(extended to 10kb at each side) would cover most of genome so it's very hard to get enough hits. So we simply shutter the whole genome and take all regions.  

### Step1, Pre-alignment   
#### 1.1, Trimming by cutadapt  
Tool: cutadapt v1.12  
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
tool: FastQC v0.11.5  
input: trimmed fastq file  
output: fastqc results  
commands: 
> fastqc -t $threads -o .  $trimed_fastq  
QC to report: fastqc results  


### Step2, reads alignment and reads distribution  
tool:   
  1, bwa v0.7.12  
  2, samtools   
  3, methylQA v0.1.9  
input: trimmed fastq file  
output: aligned bam file  
command:  
> bwa mem -t $threads  $mm10_ref_genome.fa  $trimmed.fastq | samtools view -bS - | samtools sort - -O 'bam' -o  $aln.bam -T temp_aln  
> methylQA density $chrom_size  'Trimed_rm_mapq0_chrm_'$name'.bam'
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
  
#### 4.3, PBC calculation (PCR bottlenecking coefficiency)  
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
 
 
 
