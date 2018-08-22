#!/bin/bash
# pipe usage:
# user@domain: path_to_pipe/pipe.sh -g <mm10/hg38> -r <PE/SE> -o read_file1 -p read_file2 (if PE file)
# input file: sra file, fastq file, and fastq.gz file

# pipe start
###################################################################################################
# read all necessary parameters and prepare data structure
date
pipe_version="v3.1b"
host="zhanglab/atac-seq base"

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
md5=`md5sum $0 | awk '{print $1}'`

# read parameters
while getopts m:t:g:o:p:r:i:c:h  opts
do case "$opts" in
m) marker="$OPTARG";;    # default 'no_annotation'
t) threads="$OPTARG";;    # default 24
g) species="$OPTARG";;    # hg19, hg38, mm9, mm10, danRer10
o) R1="$OPTARG";;    # PE read 1, or the SE file, or the sra file
p) R2="$OPTARG";;    # PE read 2. 
r) types="$OPTARG";;    # PE or SE;
i) ifr_parameter="$OPTARG";;  
c) methylQA_cutoff="$OPTARG";; #methylQA insertion cutoff
h) echo "
Analyzing ATAC-seq data, generating QC plot and json reports.

usage:  path-to-pipe/pipe.sh  -g <hg38/hg19/mm10/mm9/danRer10> -r <PE/SE> -o <read_file1>  -p <read_file2>

Options:    -g      input species. Please notice that each docker/singularity image is designed for one species only.
            -r      reads type, PE for paired-end reads, SE for single-end reads.
            -o      input read file1.
            -p      input read file2.
            -t      threads used for the pipe, mainly involved in cutadapt and bwa mem steps [24].
            -i      insertion free region finding parameters used by Wellington Algorithm (Jason Piper etc. 2013), see documentation for more details.
                    If you don NOT want to run IFR finding step, please just ignore the -i option; however IFR finding will use default parameters only if -i specified as 0:
                        min_lfp=5
                        max_lfp=15
                        step_lfp=2
                        min_lsh=50
                        max_lsh=200
                        step_lsh=20
                        method=BH
                        p_cutoff=0.05
                    If you want to specify your own parameter, please make sure they are in the same order and seperated by comma
                    Example: -i 5,15,2,50,200,20,BH,0.05
                    You can check the pipe log file for the parameters used by IFR code
"
exit;;
[?]) echo "
Analyzing ATAC-seq data, generating QC plot and json reports.

usage:  path-to-pipe/pipe.sh  -g <hg38/hg19/mm10/mm9/danRer10> -r <PE/SE> -o <read_file1>  -p <read_file2>

Options:    -g      input species. Please notice that each docker/singularity image is designed for one species only.
            -r      reads type, PE for paired-end reads, SE for single-end reads.
            -o      input read file1.
            -p      input read file2.
            -t      threads used for the pipe, mainly involved in cutadapt and bwa mem steps [24].
            -i      insertion free region finding parameters used by Wellington Algorithm (Jason Piper etc. 2013), see documentation for more details.
                    If you don NOT want to run IFR finding step, please just ignore the -i option; however IFR finding will use default parameters only if -i specified as 0:
                        min_lfp=5
                        max_lfp=15
                        step_lfp=2
                        min_lsh=50
                        max_lsh=200
                        step_lsh=20
                        method=BH
                        p_cutoff=0.05
                    If you want to specify your own parameter, please make sure they are in the same order and seperated by comma
                    Example: -i 5,15,2,50,200,20,BH,0.05
                    You can check the pipe log file for the parameters used by IFR code
"
exit;;
esac
done

if [ -z "$threads" ]
    then
    threads=24
fi

if [ -z "$marker" ]
    then
    marker='no_annotation'
fi

if [ -z "$methylQA_cutoff" ]
    then
    methylQA_cutoff=50
fi

if [[ $R1 == *.sra ]]
    then name=`echo ${R1%.sra}`
    echo "this is sra file, $fastq_dump_tool would be used......"
    $fastq_dump_tool $R1 --split-3
    raw1=$name'_1.fastq'
    raw2=$name'_2.fastq'
    elif [[ $R1 == *.fastq* ]] && [[ $types == PE  ]]
        then
        name=`echo ${R1%.fastq*}`
        raw1=$R1
        raw2=$R2
    elif [[ $R1 == *.fastq* ]] && [[ $types == SE ]]
        then
        name=`echo ${R1%.fastq*}`
        raw1=$R1
    else
        echo "please use fastq(or fastq.gz) file or sra file......"
        exit
fi

# analysis code
# each '0' step means those are prerequest for the following one 
# each step would assume previous steps have been processed
###################################################################################################
# step0, preparation
s0_atac_pre () {
    mkdir 'Processed_'$name
    ln -rs $R1    ./'Processed_'$name/$R1
    ln -rs $raw1  ./'Processed_'$name/$raw1  2> /dev/null
    ln -rs $raw2  ./'Processed_'$name/$raw2  2> /dev/null
    cd ./'Processed_'$name/
    source $pipe_path'/qc_source.sh' $species
    mkdir 'QC_ATAC_data_collection_'$name
    touch QC_pipe_processing.log

    # start record
    date >> QC_pipe_processing.log
    echo "Target file is $R1 $R2" >> QC_pipe_processing.log
    echo "Specified species is $species" >> QC_pipe_processing.log
    echo "types of reads is $types" >> QC_pipe_processing.log
    echo " " >> QC_pipe_processing.log
}

# Step 1.1, Trim ATAC-seq adapters and QC on seq file
s1.1_cutadapt () {
    if [[ $types == PE ]];
        then
        echo 'trimming ATAC PE reads by cutadapt'
        $cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=36  -o 'step1.1_trimed_'$name'_1.fastq' -p  'step1.1_trimed_'$name'_2.fastq'  $raw1 $raw2  > 'step1.1_'$name'_cutadapt_PE.trimlog'  
        temp=`grep "Total read pairs processed:" step1.1_*trimlog | awk '{print $5}'`  
        raw_reads=`echo ${temp//,}`  
        temp2=`grep "Pairs written" step1.1_*trimlog | awk '{print $5}'`
        written_reads=`echo ${temp2//,}`
    elif [[ $types == SE ]];
        then
        echo 'trimming ATAC SE reads by cutadapt'
        $cutadapt -a $adapter_1 --quality-cutoff=15,10 --minimum-length=36  -o  'step1.1_trimed_'$name'.fastq' $raw1 > 'step1.1_'$name'_cutadapt_SE.trimlog'  
        temp=`grep "Total reads processed:" step1.1_*trimlog | awk '{print $4}'`
        raw_reads=`echo ${temp//,}` 
        temp2=`grep "Reads written" step1.1_*trimlog | awk '{print $5}'`
        written_reads=`echo ${temp2//,}`
    fi

    if [ $? == 0 ] 
        then
        echo "step1.1, cutadapt trimming done" >> QC_pipe_processing.log
    else 
        echo "step1.1, cutadapt trimming fail......" >> QC_pipe_processing.log
        # exit 1
    fi
}

# step1.2, fastqc
s1.2_fastqc () {
    echo 'fastqc is processing fastq file......'
    [ -f ` ls 'step1.1_trimed_'$name*'.fastq' | head -1` ] && $fastqc -t $threads 'step1.1_trimed_'$name*'.fastq' -o . 
    if [ $? == 0 ] 
        then
        echo "step1.2, fastqc process done" >> QC_pipe_processing.log
    else 
        echo "step1.2, fastqc process fail......" >> QC_pipe_processing.log
        exit 1
    fi

    for zip in `ls | grep fastqc.zip`
    do
    unzip -o $zip    
    mv $zip 'step1.2_'$zip
    done

    # 1.3 fastqc data collection 
    echo -e "filename\tdeduplication_percentage\tmarker" > 'step1.3_dedup_percentage_'$name'.result'
    for file in `ls -d *fastqc/`
    do
        cd $file
        temp=`echo ${file##step1.1_trimed_}`
        out_name=`echo ${temp%*_fastqc/}`
        out_value=`grep 'Total Deduplicated Percentage' fastqc_data.txt | awk '{print $4}'`
        echo -e "$out_name\t$out_value\t$marker" >> ../'step1.3_dedup_percentage_'$name'.result'
        echo -e "item\t$out_name\t$out_name" > 'step1.3_duplication_summary_'$out_name'.result'
        grep 'Sequence Duplication Levels' -A 15 fastqc_data.txt >> 'step1.3_duplication_summary_'$out_name'.result'
        mv 'step1.3_duplication_summary_'$out_name'.result' ../'QC_ATAC_data_collection_'$name
        echo -e "$out_name\tfastqc_test" > 'step1.3_fastqc_summary_'$out_name'.result'
        awk -F "\t" '{print $1,$2}' OFS='\t' summary.txt >> 'step1.3_fastqc_summary_'$out_name'.result'
        mv 'step1.3_fastqc_summary_'$out_name'.result' ../'QC_ATAC_data_collection_'$name
        cd ..
    done

    if [ $? == 0 ] 
        then
        echo "step1.3, fastqc data_collection process done" >> QC_pipe_processing.log
    else 
        echo "step1.3, fastqc data_collection process fail......" >> QC_pipe_processing.log
    fi

    sed 1d step1.3_dedup_percentage_$name'.result' | cut -f 2 > temp_dedup.txt \
        && before_dedup=$(python -c "print(`awk '{s+=$1}END{print s}' temp_dedup.txt` * 0.01 /`cat temp_dedup.txt | wc -l`)") \
        && before_dup=$(python -c "print(1-$before_dedup*1.0)") \
        && rm temp_dedup.txt
    mv 'step1.3_dedup_percentage_'$name'.result'  ./'QC_ATAC_data_collection_'$name
    mv *fastqc* ./'QC_ATAC_data_collection_'$name

    # 1.4, get PE data R1 R2 deduplication difference percentage 
    if [[ $types == PE ]];
    then
    per1=`tail -n 2 ./'QC_ATAC_data_collection_'$name/'step1.3_dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '1p'`
    per2=`tail -n 2 ./'QC_ATAC_data_collection_'$name/'step1.3_dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '2p'`
    dif=`echo "scale=2; ($per1-$per2)*200/($per1+$per2)" | bc -l`
    else
    dif=0
    fi

    if [ $? == 0 ] 
        then
        echo "step1.4, calculate replicate difference process done" >> QC_pipe_processing.log
    else 
        echo "step1.4, calculate replicate difference process fail......" >> QC_pipe_processing.log
    fi
}

# step2.0, files check
s2.0_ref () {
    # refine chrom_size file (remove random and Unknown record)
    awk  '{ if ((length($1) < 6) && (length($1) > 1))  print $0}' OFS='\t' $chrom_size  > refined_chrom_size.txt
    chrom_size=`pwd`"/refined_chrom_size.txt"    
}

# step2.1, BWA MEM alignment
s2.1_bwa () {
    echo 'alignment by bwa......'
    $bwa mem -t $threads $bwa_ref 'step1.1_trimed_'$name*'.fastq' | $samtools view -bS - | $samtools sort - -O 'bam' -o 'step2.1_trimed_'$name'.bam' -T temp_aln
    if [ $? == 0 ] 
        then
        echo "step2.1, bwa alignment process done" >> QC_pipe_processing.log
        rm 'step1.1_trimed_'$name*'.fastq'
    else 
        echo "step2.1, bwa alignment process fail......" >> QC_pipe_processing.log
        exit 1
    fi
}

# step2.2, removing low mapQ reads and count reads distribution (mapQ=0 makes no sense here, because they are not reliable)
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
    awk -F "\t"  '{print $2}' 'step2.2_chrom_count_unique_'$name'.txt'  |  paste 'step2.2_chrom_count_'$name'.txt'  - | awk -F "\t" -v marker=$marker '{print $1,$2+0,$3+0,marker}' OFS="\t"  > ./'QC_ATAC_data_collection_'$name/'step2.2_chrom_count_'$name'.result'

    if [ $? == 0 ] 
        then
        echo "step2.2, count reads distribution process done" >> QC_pipe_processing.log
    else 
        echo "step2.2, count reads distribution process fail......" >> QC_pipe_processing.log
    fi

    rm step2.2_chrom_count*txt
}

# step2.3, preseq
s2.3_preseq () {
    $preseq lc_extrap -o 'step2.3_yield_'$name'.result' -B  'step2.1_trimed_'$name'.bam'
    if [ $? == 0 ] 
        then
        echo "step2.3, preseq lc_extrap estimate process done" >> QC_pipe_processing.log
    else 
        echo "step2.3, preseq lc_extrap estimate process fail......" >> QC_pipe_processing.log
    fi
    mv 'step2.3_yield_'$name'.result'   ./'QC_ATAC_data_collection_'$name    
}

# 3.1, methylQA
s3.1_methylQA () {
    echo 'methylQA processing......'
    echo "methylQA min insertion length choice $methylQA_cutoff" >> QC_pipe_processing.log
    grep -v ^chrM $chrom_size > nochrM_chrom_size.txt
    chrom_size=nochrM_chrom_size.txt
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
    unique_chrM_ratio=`echo "scale=4; $unique_chrM / $map_uniq" | bc -l`
    echo -e "unique_mapped\tchrM\tunique_chrM_ratio" > 'step2.2_unique_chrM_ratio_'$name'.result'
    echo -e "$map_uniq\t$unique_chrM\t$unique_chrM_ratio" >> 'step2.2_unique_chrM_ratio_'$name'.result'
    mv 'step2.2_unique_chrM_ratio_'$name'.result' ./'QC_ATAC_data_collection_'$name

    unique_no_chrM=`python -c "print($map_uniq-$unique_chrM)"`
    effect_no_chrM=`python -c "print($map_effect-$effect_chrM)"`
    nodup_ratio=`echo "scale=3; $effect_no_chrM/$unique_no_chrM" | bc -l`
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

# 3.2, normalization bedGraph -> 10M
s3.2_nomral_bg () {
    echo 'normalization bedGraph......'

    # add a new bigwig file without black list
    intersectBed -iobuf 200M -a 'step3.1_methylQA_'*$name'.open.bedGraph'  -b $black_list -v  > rmbl.bedGraph
    bedGraphToBigWig rmbl.bedGraph $chrom_size 'step3.2_rmbl_'$name'.bigWig' && rm 'step3.1_methylQA_'*$name'.open.bedGraph'

    # normalization
    norm=`grep 'non-redundant'  step3.1_methylQA_*report | awk '{print $6}'`
    factor=`echo "scale=3; $norm/10000000" | bc -l`
    awk -v factor=$factor '{print $1,$2,$3,$4/factor}' OFS='\t' rmbl.bedGraph  >  'step3.2_normalized_per_10M_'$name'.open.bedGraph'
    bedGraphToBigWig  'step3.2_normalized_per_10M_'$name'.open.bedGraph'   $chrom_size  'step3.2_normalized_per_10M_'$name'.bigWig' && rm 'step3.2_normalized_per_10M_'$name'.open.bedGraph'
    rm rmbl.bedGraph

    if [ $? == 0 ] 
        then
        echo "step3.2, normalization process done" >> QC_pipe_processing.log
    else 
        echo "step3.2, normalization process fail......" >> QC_pipe_processing.log
    fi
}

# 3.3, peak calling
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

# step4.0, set variable
s4.0_set () {
    peak='step3.4_peakcall_'$name'_peaks.narrowPeak'
    bed='step3.3_rmbl_'$name'.open.bed'    
}

# 4.1, RUP and insertion site
s4.1_rup () {
    total=`wc -l $bed |awk '{print $1}'`
    sum=`intersectBed -iobuf 200M -a $bed -b $peak -f 0.5 -u | wc -l`
    ratio=`echo "scale=2; $sum*100/$total" | bc -l`
    if [ $? == 0 ] 
        then
        echo "step4.1, reads unpder peak ratio calculation process done" >> QC_pipe_processing.log
    else 
        echo "step4.1, reads unpder peak ratio calculation process fail......" >> QC_pipe_processing.log
    fi

    # 4.1.2, add insertion site bigwig
    awk '{mid=int(($3+$2)/2); if($6=="+") {print $1"\t"mid"\t"mid+1"\t"1} else {print $1"\t"mid-1"\t"mid"\t"1}}'  \
     $bed |  sort -k1,1 -k2,2n | uniq -c | awk -F " " '{print $2"\t"$3"\t"$4"\t"$1}' > step4.2_insertion_site_$name.bedGraph
    bedGraphToBigWig  step4.2_insertion_site_$name.bedGraph  $chrom_size  step4.2_insertion_site_$name'.bigWig'
    if [ $? == 0 ] 
        then
        echo "step4.1.2, insertion site process done" >> QC_pipe_processing.log
    else 
        echo "step4.1.2, insertion site process fail......" >> QC_pipe_processing.log
    fi    
}

# 4.2, enrichment
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


# 4.4, saturation analysis
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

# 4.5, background
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

# step 4.6, visualization
s4.6_visualization () {
    #summarize results
    echo -e "file\ttotal\twritten_reads\tmapped\tmapped_ratio\tuniq_mapped\tnon_redundant_uniq_mapped\teffect_ratio\tfastqc_dup\tafter_align_dup\tnumber_of_reads_under_peak\trup_ratio\tsub10M_enrichment\tcoding_enrichment\tbg_gt37_percentage" >  'QC_data_collection_'$name'.result'
    echo -e "$name\t$raw_reads\t$written_reads\t$map_mapped\t$mapped_ratio\t$map_uniq\t$map_effect\t$effect_ratio\t$before_dup\t$after_dup\t$sum\t$ratio\t$enrichment\t$enrichment_ratio\t$ra_more" >>  'QC_data_collection_'$name'.result'
    mv 'QC_data_collection_'$name'.result'  ./'QC_ATAC_data_collection_'$name

    if [ -z $name ] || [ -z $raw_reads ] || [ -z $map_mapped ] || [ -z $mapped_ratio ] || [ -z $map_uniq ] || [ -z $written_reads ] || [ -z $map_effect ] || [ -z $effect_ratio ] || [ -z $nodup_ratio ] || [ -z $sum ]|| [ -z $ratio ]|| [ -z $dif ] || [ -z $before_dedup ]
    then
        echo "step4.6, sumarizing result process fail......" >> QC_pipe_processing.log
    else
        echo "step4.6, sumarizing result process done" >> QC_pipe_processing.log
    fi   

    # plot and json
    time=`head -1 QC_pipe_processing.log | sed 's/ /_/g'`
    image_id=`bash $pipe_path'/find_image_ID_digest.sh' $host  2> /dev/null | awk '{print $2}'`
    if [ -z "$image_id" ]
    then
    image_id="failed_to_get_id"
    fi

    # clean result
    find . -name "*.result" | xargs sed -i 's/^-e //'
    cd ./'QC_ATAC_data_collection_'$name
    Rscript $pipe_path'/visualization.R' $name $pipe_path'/../atac_ref/mm10_encode_pe'  $species  $written_reads $unique_chrM_ratio $pipe_version $time $image_id 
    if [ $? == 0 ] 
        then
        echo "step4.6, plot process done" >> ../QC_pipe_processing.log
    else 
        echo "step4.6, plot process fail......" >> ../QC_pipe_processing.log
    fi

    sed 's/\[/{/g' $name'_report.json' | sed '/      {/d' | sed '/\]/d' |\
        sed 's/      }/    },/g' | sed 's/"!/{/g' | sed 's/!"/}/g' | sed 's/"?/[/g' | sed 's/?"/]/g' |\
        sed 's/@/"/g' | tac | sed '3s/},/}/g' | sed '1,2d' | tac | cat - <(echo "  },") <(sed '1d' $pipe_path'/../atac_ref/mm10_encode_pe/encode_pe.json') | sed 's/\\r//g' | sed "s/MD5ToBeChange/$md5/g" > QC_$name'.json'
    rm $name'_report.json'
    mv QC_$name'.json' ../
    paste <(cut -f 1-8 QC_data_collection_${name}.result) <(cut -f 5 step3.1_useful_reads_${name}.result) <(cut -f 9- QC_data_collection_${name}.result) > QC_table_${name}.result

    mkdir 'plots_collection_'$name
    mv *png 'plots_collection_'$name
    rm $name'_report.txt' 
    cd ..

    rm pesudo_bl.txt 2> /dev/null
    rm refined_chrom_size.txt 
    find -type l -delete

    multiqc .
    rename 's/multiqc/step4.6_multiqc/' multiqc*
}

s4.7_ifr_finding () {
    # set default parameters:
    if (( $ifr_parameter == 0 )); then
        min_lfp=5
        max_lfp=15
        step_lfp=2
        min_lsh=50
        max_lsh=200
        step_lsh=20
        method=BH
        p_cutoff=0.05
    else
        min_lfp=`cut -d"," -f 1 <<< $ifr_parameter`
        max_lfp=`cut -d"," -f 2 <<< $ifr_parameter`
        step_lfp=`cut -d"," -f 3 <<< $ifr_parameter`
        min_lsh=`cut -d"," -f 4 <<< $ifr_parameter`
        max_lsh=`cut -d"," -f 5 <<< $ifr_parameter`
        step_lsh=`cut -d"," -f 6 <<< $ifr_parameter`
        method=`cut -d"," -f 7 <<< $ifr_parameter`
        p_cutoff=`cut -d"," -f 8 <<< $ifr_parameter`
    fi

    echo "step4.7, IFR finding parameters used:
    min_lfp is $min_lfp
    max_lfp is $max_lfp
    step_lfp is $step_lfp
    min_lsh is $min_lsh
    max_lsh is $max_lsh
    step_lsh is $step_lsh
    method is $method
    p_cutoff is $p_cutoff
    "  >> QC_pipe_processing.log

    awk '{print $1"\t"$2-50"\t"$3+50"\t"$4}' "step3.4_peakcall_"$name"_peaks.narrowPeak" |\
    awk '{if(($3-$2)>1000) {printf $1"\t"$2"\t"; printf "%.0f",($2+$3)/2; print "\t"$4".1"; printf $1"\t"; printf "%.0f",($2+$3)/2; print "\t"$3"\t"$4".2"} else print $0}' |\
    awk '{if(($3-$2)>1000) {printf $1"\t"$2"\t"; printf "%.0f",($2+$3)/2; print "\t"$4".1"; printf $1"\t"; printf "%.0f",($2+$3)/2; print "\t"$3"\t"$4".2"} else print $0}' > temp.peak

    split -n l/$threads temp.peak
    rm temp.peak

    [ -f temp.txt ] && rm list
    for file in `ls xa*`
    do
        intersectBed -a "step4.2_insertion_site_"$name".bedGraph" -b $file -wa -wb | uniq > $file".bed"
        echo "Rscript $pipe_path'/ATAC-seq_wellington.R' "$file".bed IFR_"$file".txt" $min_lfp $max_lfp $step_lfp $min_lsh $max_lsh $step_lsh $method $p_cutoff >> list
        rm $file
    done

    cat list | parallel 
    if [ $? == 0 ] 
        then
        echo "step4.7, IFR finding process done" >> QC_pipe_processing.log
    else 
        echo "step4.7, cIFR finding process fail......" >> QC_pipe_processing.log
    fi

    rm xa*bed list
    cat IFR*txt | sed "s/\"//g" | sort -k1,1V -k2,2n | awk '{print $1"\t"$2"\t"$3"\t""found_IFR_"NR"\t"$4"\t"".""\t"$5}' > "step4.7_IFR_"$name".bed"
    rm IFR*txt
}

# run pipe ()
###################################################################################################
# step-by-step
s0_atac_pre
s1.1_cutadapt
s1.2_fastqc
s2.0_ref
s2.1_bwa
s2.2_distri
s2.3_preseq  
s3.1_methylQA
s3.2_nomral_bg
s3.3_peakcall
s4.0_set
s4.1_rup 
s4.2_enrich
s4.4_saturation 
s4.5_background  
s4.6_visualization 

if [ -z "$ifr_parameter" ]; then
    echo "step4.7 ifr finding is ommited" >> QC_pipe_processing.log
else
    s4.7_ifr_finding
fi

echo "Processing $name done"
echo "Processing $name done"
echo "Processing $name done"
cd ..
date
















