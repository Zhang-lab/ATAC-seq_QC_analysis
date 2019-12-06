# Documentation v1.00
ATAC-seq quality control and analysis pipeline for Bo Zhang's lab  
Last edit: 02/11/2019  
For any question please contact: shaopeng.liu@wustl.edu  

**Outline**  
I, Output example and annotation  
II, Downstream analysis  
III, Visulization through qATACViewer  
IV, Terms  
V, Data processing details  

&nbsp;
&nbsp;
&nbsp;
&nbsp;
## I, Output example  
After running the pipeline, there will be a folder called **Processed_${name}**, all intermediate files and final output files are stored there. And it looks like this (using an ENCODE data for example):  
![output image](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/ATAC_out.png)  

&nbsp;
&nbsp;
&nbsp;
&nbsp;

There should be 15 files in total, and the contents are listed below.  
 
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

&nbsp;
&nbsp;
&nbsp;
&nbsp;

If you happens to open the first folder which record intermediate files, you will see:  
![QC image](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/QC_col_folder_example.png)  

There will be 33 files, but don't worry, all of them are summaized together in the json file. You might be interested some figures stored in the first folder "plots_collection_${name}". If you open it, you will get this.
![plot collection](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/plot_example.png)  

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


&nbsp;
&nbsp;
&nbsp;
&nbsp;
## II, Downstream analysis  
1. merge QC data for a batch of files  
```
# come to the parent folder where you store all those files

singularity exec <path-2-singularity-image> bash /atac_seq/pipe_code/batch_collection_v3.1.sh
```


2. Differential Accessible Region (DAR) detection
Borrowing the idea of Differentially Expressed Genes (DEG) analysis, we utilize **DESeq2** to identify peak regions in ATAC-seq with differentially insertion events. 
```
# similar to DEG, you will need 2 conditions and sample size for each condition is 2, unbalanced design is supported.
# please put bed and peak file (*step3.3_rmbl_${name}.open.bed*, and *step3.4_peakcall_${name}_peaks.narrowPeak*) of each data into same folder
# please make sure for each condition, there is a keyword to group them together, e.g. group1 group2 as prefix. You can add it if necessary.

singularity exec <path-2-singularity-image> bash /atac_seq/pipe_code/DOR_analysis.sh ${group1_key} ${group2_key}
```


&nbsp;
&nbsp;
&nbsp;
&nbsp;
## III, Visulization through qATACViewer
To make better visulization, we have also prepared a tool named **qATACViewer**, please **[click here](https://github.com/lidaof/qATACviewer/tree/localjson)** to find its github page.  

And please **[click here](http://brc.wustl.edu/SPACE/shaopengliu/atac_v1/20190211_AIAP_documentation_demo/QAViewer_ENCODE_liver.pdf)** for examples.

**Usage**
1. copy the git repository to local: `git clone -b localjson https://github.com/lidaof/qATACviewer.git` (please makesure to use the `localjson` branch)  
2. update the path in `./qATACViewer/frontend/src/data.json` **(we are improving this part)**    
3. follow the instruction on the the **Readme** of it  

&nbsp;
&nbsp;
&nbsp;
&nbsp;
## IV, Term and definition  
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

&nbsp;
&nbsp;
&nbsp;
&nbsp;
## V, Data Processing details  
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
```
$cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=36  -o 'step1.1_trimed_'$name'_1.fastq' -p  'step1.1_trimed_'$name'_2.fastq'  $input1 $input2  > 'step1.1_'$name'_cutadapt_PE.trimlog'  
```
  
#### 1.2, Fastq file quality assessment  
tool: FastQC v0.11.7  
input: trimmed fastq file  
output: fastqc results  
```
[ -f ` ls 'step1.1_trimed_'$name*'.fastq' | head -1` ] && $fastqc -t $threads 'step1.1_trimed_'$name*'.fastq' -o . 
```


### Step2, reads alignment and reads distribution  
#### 2.1, BWA MEM aligner  
tool: bwa v0.7.12  
input: trimed fastq file  
```
$bwa mem -t $threads $bwa_ref 'step1.1_trimed_'$name*'.fastq' | $samtools view -bS - | $samtools sort - -O 'bam' -o 'step2.1_trimed_'$name'.bam' -T temp_aln
```

#### 2.2, reads distribution  
tool: samtools 1.3.1, methylQA v0.2.1    
input: aligned bam file  
output: reads distribution count  
```
s2.2_distri () {
    $samtools view -h 'step2.1_trimed_'$name'.bam' > input.sam \
        && awk '$5>0' input.sam | sed '/^@/d' - | cat <(grep '^@' input.sam) - > output.sam \
        && cat output.sam | awk '{print $3}' |  sort  -k1,1V |   uniq -c > count_no_mapq0.txt \
        && awk '! /random/ && ! /Un/ && /chr/  ' count_no_mapq0.txt  | awk '{print $2, $1}'  OFS="\t"  | sort  -k1,1 -V -s > temp2.2.txt 

    if [[ $types == SE ]]; then
    mv temp2.2.txt  'step2.2_chrom_count_'$name'.txt'
    elif [[ $types == PE ]]; then
    awk '{$2=int($2*0.5); print}' OFS="\t" temp2.2.txt > 'step2.2_chrom_count_'$name'.txt' && rm temp2.2.txt
    fi

    # only effect reads
    $methylQA density -S $chrom_size  output.sam
    cut -f 1 output.extended.bed | uniq -c >  count_unique.txt
    awk '{print $2, $1}'  OFS="\t" count_unique.txt | sort  -k1,1 -V -s > 'step2.2_chrom_count_unique_'$name'.txt'
    effect_chrM=`grep chrM output.extended.bed | wc -l`

    # get chrM count in uniquely mapped reads
    # to keep the count consistent, the results are all from methylQA
    # when count directly from output.sam file, please pay attention to the unpaired reads
    cat <(samtools view -H input.sam) <(awk '$3=="chrM"' input.sam)  | $methylQA density -S -r -o temp $chrom_size  -
    unique_chrM=`grep chrM temp.extended.bed | wc -l`

    rm temp*
    rm output*
    rm count*.txt
    rm input.sam
    awk -F "\t"  '{print $2}' 'step2.2_chrom_count_unique_'$name'.txt'  |  paste 'step2.2_chrom_count_'$name'.txt'  - | awk -F "\t" -v marker=$marker '{print $1,$2+0,$3+0,marker}' OFS="\t" | awk 'length($1)<7'  > ./'QC_ATAC_data_collection_'$name/'step2.2_chrom_count_'$name'.result'

    if [ $? == 0 ] 
        then
        echo "step2.2, count reads distribution process done" >> QC_pipe_processing.log
    else 
        echo "step2.2, count reads distribution process fail......" >> QC_pipe_processing.log
    fi

    rm step2.2_chrom_count*txt
}
``` 

#### 2.3, library complexity estimate  
tool: preseq 2.0.0 (from subread)  
input: aligned bam file  
```
$preseq lc_extrap -o 'step2.3_yield_'$name'.result' -B  'step2.1_trimed_'$name'.bam'
```

### Step3, reads filter and peak call  
tool:   
  1, methylQA v0.1.9  
  2, macs2 v2.1  
input:   
  1, alined file as $aln.bam  
  2, effect reads file from 1 as $reads.open.bed  
output: peak file  
```
s3.1_methylQA () {
    echo 'methylQA processing......'
    echo "methylQA min insertion length choice $methylQA_cutoff" >> QC_pipe_processing.log
    grep -v ^chrM $chrom_size > nochrM_chrom_size.txt \
        && rm $chrom_size \
        && chrom_size=nochrM_chrom_size.txt
    echo "the chrom file for methylQA is $chrom_size" >> QC_pipe_processing.log
    $methylQA atac -X $methylQA_cutoff  -o step3.1_methylQA_$name  $chrom_size  'step2.1_trimed_'$name'.bam'

    if [ $? == 0 ] 
        then
        echo "step3.1, mathylQA atac process done" >> QC_pipe_processing.log
    else 
        echo "step3.1, mathylQA atac process fail......" >> QC_pipe_processing.log
    fi

    # mapping status
    map_mapped=`grep 'mappable reads' step3.1_methylQA_$name'.report' | awk '{print $4}'`
    map_uniq=`grep '(mapQ >= 10)' step3.1_methylQA_$name'.report' | awk '{print $8}'`
    map_effect=`grep 'non-redundant'  step3.1_methylQA_$name'.report' | awk '{print $6}'`
    mapped_ratio=`echo "scale=2; $map_mapped/$raw_reads" | bc -l`
    effect_ratio=`echo "scale=2; $map_effect/$raw_reads" | bc -l`

    # unique chrM ratio from step2.2
    unique_chrM_ratio=`python -c "print($unique_chrM*1.0 / ($unique_chrM+$map_uniq) )"`
    echo -e "non_chrM_unique_mapped\tchrM\tunique_chrM_ratio" > 'step2.2_unique_chrM_ratio_'$name'.result'
    echo -e "$map_uniq\t$unique_chrM\t$unique_chrM_ratio" >> 'step2.2_unique_chrM_ratio_'$name'.result'
    mv 'step2.2_unique_chrM_ratio_'$name'.result' ./'QC_ATAC_data_collection_'$name

    #unique_no_chrM=`python -c "print($map_uniq-$unique_chrM)"`
    #effect_no_chrM=`python -c "print($map_effect-$effect_chrM)"`
    nodup_ratio=`echo "scale=3; $map_effect/$map_uniq" | bc -l`
    after_dup=$(python -c "print(1-$nodup_ratio*1.0)")

    useful=`grep 'non-redundant'  step3.1_methylQA_*.report | awk '{print $6}'`
    single_end=`wc -l *open.bed | awk '{print $1}'`
    uf_ratio=`echo "scale=3; $useful / $raw_reads" | bc -l`
    echo -e "file\ttotal\tuseful\tuseful_ratio\tsingle_end" > 'step3.1_useful_reads_'$name.result
    echo -e "$name\t$raw_reads\t$useful\t$uf_ratio\t$single_end" >> 'step3.1_useful_reads_'$name.result
    mv 'step3.1_useful_reads_'$name.result  ./'QC_ATAC_data_collection_'$name
    sort -n 'step3.1_methylQA_'*$name'.insertdistro' | uniq -c | awk '{print $2,$1}' > 'step3.1_insertion_distri_'$name'.result' && rm 'step3.1_methylQA_'*$name'.insertdistro'
    mv 'step3.1_insertion_distri_'$name'.result'  ./'QC_ATAC_data_collection_'$name
    rm step3.1_methylQA_*bigWig
}

s3.3_peakcall () {
    echo 'peak calling......'

    awk '{if ($2 > $3)sub($2, 0); print}' OFS="\t" 'step3.1_methylQA_'$name'.open.bed' > temp.open.bed \
    && intersectBed -iobuf 200M  -a temp.open.bed  -b $black_list   -v > 'step3.3_rmbl_'$name'.open.bed' \
    && rm temp.open.bed 'step3.1_methylQA_'$name'.open.bed'

    $macs2 callpeak -t 'step3.3_rmbl_'$name'.open.bed' -g $macs2_genome -q 0.01 -n 'step3.4_peakcall_'$name  --keep-dup 1000 --nomodel --shift 0 --extsize 150

    if [ $? == 0 ] 
        then
        echo "step3.4, macs2 peak calling process done" >> QC_pipe_processing.log
    else 
        echo "step3.4, macs2 peak calling process fail......" >> QC_pipe_processing.log
    fi

    mv step3.1_methylQA_$name*  'QC_ATAC_data_collection_'$name

    # peak length distribution:
    awk '{print $3-$2+1}' 'step3.4_peakcall_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' > 'step3.4_peak_length_distri_'$name'.result'
    mv 'step3.4_peak_length_distri_'$name'.result'  ./'QC_ATAC_data_collection_'$name    
}
```
note: methylQA atac would do a PE->SE transformation based on insertion site for PE data.  


### Step4, Post alignment analysis:  
#### 4.1, Reads under peak ratio  
tool: bash
input: reads file as $bed, and peak file as $peak  
output: reads under peak ratio  
```
total=`wc -l $bed |awk '{print $1}'`
sum=`intersectBed -iobuf 200M -a $bed -b $peak -f 0.5 -u | wc -l`
ratio=`echo "scale=2; $sum*100/$total" | bc -l`
```

#### 4.2, Calculate enrichment ration of top 20k peaks against overall background  
tool: bash  
input: reads enrichment in coding promoters, and normalized enrichment for all peak region
output: enrichment ratio file  
```
s4.2_enrich () {
    # 4.2.1, new enrichment from RUP based on 10M sub-sampling with adjustment
    # numerator =  ($rupn+10000000*$peak_length / $genome_size) / $peak_length
    # denominator = (10000000+$useful_ends) / ($genome_size-$peak_length))
    ## $1 for original open.bed, $2 for useful_single ends/sub-sample 10M
    cal_enrich () {
        shuf $1 | head -10000000  > temp.open.bed
        $macs2 callpeak -t temp.open.bed  -g $macs2_genome -q 0.01 -n temp_peak  --keep-dup 1000 --nomodel --shift 0 --extsize 150
        peak_length=`awk '{s+=$3-$2+1}END{print s}' 'temp_peak_peaks.narrowPeak'`
        rupn=`intersectBed -iobuf 200M -a temp.open.bed -b 'temp_peak_peaks.narrowPeak' -f 0.5 | wc -l`
        upper=`python -c "print(1.0*($rupn+10000000*$peak_length/$genome_size)/$peak_length)"`
        lower=`python -c "print(1.0*($2+10000000)/($genome_size-$peak_length))"`
        enrichment=`python -c "print(1.0*$upper/$lower)"`
        rup=`python -c "print(1.0*$rupn/10000000)"`
        echo -e "name\trupn\trup\tcoverage\tenrichment" > 'step4.2_sub10M_enrichment_'$name'.result'
        echo -e "$name\t$rupn\t$rup\t$peak_length\t$enrichment" >> 'step4.2_sub10M_enrichment_'$name'.result'
        mv 'step4.2_sub10M_enrichment_'$name'.result'  ./'QC_ATAC_data_collection_'$name
        rm temp.open.bed temp_peak_*
    }

    total=`wc -l $bed |awk '{print $1}'`
    if (( $total > 10000000 )) 
    then
        cal_enrich $bed 10000000
    else
        cal_enrich $bed $total
        echo "Warning: the open.bed file contains less than 10M reads" >> QC_pipe_processing.log
        echo "Warning: the enrichment is calculated by the original bed file, and may not be reliable" >> QC_pipe_processing.log
    fi

    if [ $? == 0 ] 
        then
        echo "step4.2.1, sub10M enrichment ratio process done" >> QC_pipe_processing.log
    else 
        echo "step4.2.1, sub10M enrichment ratio process fail......" >> QC_pipe_processing.log
    fi

    # 4.2.2, coding promoter enrichment
    # coding enrichment = ( reads in promoter / promoter length)  /  (total reads / genome size)
    denominator=`echo "scale=10; $total / $genome_size" | bc -l`
    intersectBed -iobuf 200M -a $peak -b $coding_promoter -u > promoter_peak.bed
    reads_in_promoter=`intersectBed -iobuf 200M -a $bed -b promoter_peak.bed -f 0.5 -u | wc -l | awk '{print $1}'`
    promoter_number=`intersectBed -iobuf 200M -a $coding_promoter -b promoter_peak.bed -F 0.5 -u | wc -l | awk '{print $1}'`
    promoter_length=`echo "$promoter_number * 2000+0.001" | bc -l`  
    enrichment_ratio=`echo "scale=3; $reads_in_promoter / $promoter_length / $denominator" | bc -l`
    if [ $? == 0 ] 
        then
        echo "step4.2.2, coding promoter enrichment ratio process done" >> QC_pipe_processing.log
    else 
        echo "step4.2.2, coding promoter enrichment ratio process fail......" >> QC_pipe_processing.log
    fi
    echo -e "name\ttotal_reads\tpromoter_number\treads_in_promoter\tenrichment_ratio" > 'step4.2_enrichment_ratio_in_promoter_'$name'.result'
    echo -e "$name\t$total\t$promoter_number\t$reads_in_promoter\t$enrichment_ratio" >> 'step4.2_enrichment_ratio_in_promoter_'$name'.result'
    mv 'step4.2_enrichment_ratio_in_promoter_'$name'.result'  'QC_ATAC_data_collection_'$name
    rm promoter_peak.bed
}
``` 
  
#### 4.4, Saturation analysis  
tool: macs2 v2.1  
input: effect reads  
commands:  
> randomly sampling reads file from 5%, 10%, 20% to 90% of reads file  
> call peak and calculate the peaks in sub-sampling file  
```
s4.4_saturation () {
    # subsampling:
    total=`wc -l $bed |awk '{print $1}'`
    for number in 5 10 20 30 40 50 60 70 80 90
    do
        sample_ratio=$(( total * $number / 100 ))
        shuf $bed | head -$sample_ratio >  'Trimed_rmbl_'$name'_sample'$number'.open.bed'
    done

    if [ $? == 0 ] 
        then
        echo "step4.4, saturation subsampling process done" >> QC_pipe_processing.log
    else 
        echo "step4.4, saturation subsampling process fail......" >> QC_pipe_processing.log
    fi

    # call peak
    mkdir saturation_$name
    mv *sample*.open.bed ./saturation_$name/
    ln -rs $bed ./saturation_$name/
    cp step3.4_peakcall_*Peak ./saturation_$name/'peakcall_Trimed_rmbl_'$name'.open.bed_peaks.narrowPeak'
    cd ./saturation_$name
    for file in `ls 'Trimed_rmbl_'$name'_sample'*'.open.bed'`;
    do
        $macs2 callpeak -t $file -g $macs2_genome -q 0.01 -n 'peakcall_'$file --keep-dup 1000 --nomodel --shift 0 --extsize 150
    done

    if [ $? == 0 ] 
        then
        echo "step4.5, saturation call peak process done" >> ../QC_pipe_processing.log
    else 
        echo "step4.5, saturation call peak process fail......" >> ../QC_pipe_processing.log
    fi

    echo "peak calling done......"

    # summarise results
    echo -e "5\n`seq 10 10 100`" > saturation_points.txt

    for file in `ls *open.bed`
    do
    read_num=`wc -l $file | awk '{print $1}'`
    echo `echo "scale=2; $read_num / 1000000" | bc -l`>> temp44.txt 
    rm $file
    done
    sort -k1,1n temp44.txt > saturation_reads.txt
    rm temp44.txt

    total_region=`awk '{s+=$3-$2+1}END{print s}' 'peakcall_Trimed_rmbl_'$name'.open.bed_peaks.narrowPeak'`


    for number in 5 10 20 30 40 50 60 70 80 90
    do
    file='peakcall_Trimed_rmbl_'$name'_sample'$number'.open.bed_peaks.narrowPeak'
    peak_number=`wc -l $file | awk '{print $1}'`
    peak_region=`intersectBed -iobuf 200M -a $file -b 'peakcall_Trimed_rmbl_'$name'.open.bed_peaks.narrowPeak' | awk '{s+=$3-$2+1}END{print s}'`
    if [ -z "$peak_region" ]; then
    peak_region=0
    fi
    echo `echo "scale=2; $peak_region / $total_region" | bc -l` >> temp443.txt
    echo $peak_number >> temp442.txt
    done

    if [ $? == 0 ] 
        then
        echo "step4.5, saturation results collection process done" >> ../QC_pipe_processing.log
    else 
        echo "step4.5, saturation results collection process fail......" >> ../QC_pipe_processing.log
    fi

    echo `wc -l 'peakcall_Trimed_rmbl_'$name'.open.bed_peaks.narrowPeak' | awk '{print $1}'` >> temp442.txt
    mv temp442.txt  saturation_peak.txt
    echo 1 >> temp443.txt
    mv temp443.txt  saturation_ratio.txt
    paste saturation_points.txt  saturation_reads.txt  saturation_peak.txt   saturation_ratio.txt  > temp444.txt

    echo -e "file\t$name'_read'\t$name'_peak'\t$name'_ratio'\tmarker"  > 'step4.4_saturation_'$name'.result'  
    awk -v marker=$marker '{print $0,marker}' OFS='\t' temp444.txt >> 'step4.4_saturation_'$name'.result'  
    rm temp444.txt
    rm saturation*.txt
    mv 'step4.4_saturation_'$name'.result'  ../'QC_ATAC_data_collection_'$name
    cd ..    
    mv saturation_$name 'QC_ATAC_data_collection_'$name
}
```
  
#### 4.5, Calculate background  
tool:   
  1, Python  
  2, intersectBed  
input: random regions from genome that are at least 10kb from any peak  
command:  
> generate random regions from genome  
> remove those with peak in 10kb range  
> calculate RPKM for those region  
```
s4.5_background () {
    # the exit singal 1 happens then peak number < 100, thus the 2 "mv bin.txt" command would have minor error, doesn't influence results
    # signal part
    intersectBed -iobuf 200M -a $peak -b $promoter_file -u | awk '{print $1"\t"$2"\t"$3"\t""1""\t"$9}' > promoter.narrowPeak
    intersectBed -iobuf 200M -a $peak -b $promoter_file -v | awk '{print $1"\t"$2"\t"$3"\t""0""\t"$9}' > non-promoter.narrowPeak

    echo -e "num_peaks_in_promoter\tnum_peaks_in_non-promoter\tnum_reads_in_promoter_peaks\tnum_reads_in_non-promoter_peaks" > 'promoter_percentage_'$name'.result'

    peak1=`wc -l promoter.narrowPeak | awk '{print $1}'`
    peak2=`wc -l non-promoter.narrowPeak | awk '{print $1}'`
    read1=`intersectBed -iobuf 200M -a $bed -b promoter.narrowPeak -u -f 0.50 | wc -l`
    read2=`intersectBed -iobuf 200M -a $bed -b non-promoter.narrowPeak -u -f 0.50 | wc -l`

    echo -e "$peak1\t$peak2\t$read1\t$read2" >> 'promoter_percentage_'$name'.result'
    sed -i 's/^-e //' 'promoter_percentage_'$name'.result'

    cat promoter.narrowPeak non-promoter.narrowPeak | sort -k5 -n -r > top10k.narrowPeak

    if (( `cat top10k.narrowPeak | wc -l` > 100 ))
    then
    python $pipe_path'/promoter_bin.py' top10k.narrowPeak
    else
    echo "Warning: total peak is fewer than 100, promoter bin step would be skipped. At least 100 peaks are required." >> ./QC_pipe_processing.log
    fi

    rm promoter.narrowPeak
    rm non-promoter.narrowPeak
    rm top10k.narrowPeak

    # background noise part
    awk  '{ if ((length($1) < 6) && (length($1) > 1))  print $0}' OFS='\t' $chrom_size > temp45.txt
    python2.7  $pipe_path'/random_chr.py' temp45.txt
    size=`wc -l $bed | awk '{print $1}'`

    awk '{print $1"\t"int(($3+$2)/2)-100000"\t"int(($3+$2)/2)+100000"\t"$4}' $peak > temp45
    awk '{if ($2<0) $2=0; print $0}' OFS="\t" temp45 > temp452
    mv temp452 temp45
    intersectBed -iobuf 200M -a chr.peak -b temp45 -v | shuf - | head -50000 | sort -k1,1V -k2,2n > background
    intersectBed -iobuf 200M -a $bed -b background -u -f 0.5 | sort -k1,1V -k2,2n > temp45
    python $pipe_path'/rpkm_bin.py' background temp45 $size

    if [ $? == 0 ] 
        then
        echo "step4.6, background evaluation process done" >> ./QC_pipe_processing.log
    else 
        echo "step4.6, background evaluation process fail......" >> ./QC_pipe_processing.log
    fi

    mv reads.txt 'background_'$name'.result'
    rm temp45
    rm background
    mv bin.txt 'bin_'$name'.result' 2> /dev/null

    bg_total=`wc -l background*.result | awk '{print $1}'`  
    bg_half_thres=`awk '$6<=0.188 {print $0}' background*.result | wc -l`  
    bg_less=`awk '$6<=0.377 {print $0}' background*.result | wc -l`  
    bg_more=`awk '$6>0.377 {print $0}' background*.result | wc -l` 
    ra_half_thres=`echo "scale=2; $bg_half_thres*100 / $bg_total" | bc -l`  
    ra_less=`echo "scale=2; $bg_less*100 / $bg_total" | bc -l`  
    ra_more=`echo "scale=2; $bg_more*100 / $bg_total" | bc -l`  
    echo -e "$ra_half_thres\t$ra_less\t$ra_more" > 'dichoto_bg_'$name'.result'  

    mv 'dichoto_bg_'$name'.result'   ./'QC_ATAC_data_collection_'$name/'step4.5_dichoto_bg_'$name'.result'
    mv background_$name'.result'  ./'QC_ATAC_data_collection_'$name/step4.5_background_$name'.result'
    mv promoter_percentage_$name'.result' ./'QC_ATAC_data_collection_'$name/step4.5_promoter_percentage_$name'.result'
    mv bin_$name'.result' ./'QC_ATAC_data_collection_'$name/step4.5_bin_$name'.result'  2> /dev/null
    rm temp45.txt chr.peak
}
```
   
#### 4.6, plot all results with reference dataset collection (currently ENCODE mm10 PE data)  
tool: R  
QC to report: all plots and report file 
**[click here](https://github.com/Zhang-lab/ATAC-seq_QC_analysis/blob/master/pipe_code/visualization.R)** for code
 
 
