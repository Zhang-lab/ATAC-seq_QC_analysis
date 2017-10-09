#!/bin/bash
# pipe usage:
# user@domain: path_to_pipe/pipe.sh -g -r <PE/SE> -o read_file1 -p read_file2 (if PE file)
# Optional parameter: -t -m -h for help
# -t for threads number, default 24
# -m for marker, default 'unmarked'
# -h for help
# input file: sra file, fastq file, and fastq.gz file

# pipe start
###################################################################################################
# Preparation:
date

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 

# read parameters
while getopts m:t:g:o:p:r:  opts
do case "$opts" in
m) marker="$OPTARG";;	# default 'unmarked'
t) threads="$OPTARG";;	# default 24
g) species="$OPTARG";;	# hg19, hg38, mm9, mm10, danRer10
o) R1="$OPTARG";;    # PE read 1, or the SE file, or the sra file
p) R2="$OPTARG";;    # PE read 2. 
r) types="$OPTARG";;	# PE or SE;
h) echo "usage:  path-to-pipe/pipe.sh  -g <hg38/hg19/mm10/mm9/danRer10/personal>  -r <PE/SE> -o read_file1  -p read_file2 (if necessary)"
exit;;
[?]) "Usage: ./pipe.sh  -g hg38/hg19/mm10/mm9/danRer10 -o read_file1  -r PE/SE";;
esac
done

if [ -z "$threads" ]
then
threads=24
fi

if [ -z "$marker" ]
then
echo "you didn't specify the marker, are you sure to keep the default unmarked"
marker='unmarked'
fi

source $pipe_path'/qc_source.sh' $species

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

mkdir 'Processed_'$name
mv $R1    ./'Processed_'$name/
mv pesudo_bl.txt  ./'Processed_'$name/  2> /dev/null
mv $raw1  ./'Processed_'$name/  2> /dev/null
mv $raw2  ./'Processed_'$name/  2> /dev/null
cd ./'Processed_'$name/
mkdir 'data_collection_'$name
touch pipe_processing.log

# refine chrom_size file (remove random and Unknown record)
awk  '{ if ((length($1) < 6) && (length($1) > 1))  print $0}' OFS='\t' $chrom_size  > refined_chrom_size.txt
chrom_size=`pwd`"/refined_chrom_size.txt"

# start record
date >> pipe_processing.log
echo "Target file is $R1 $R2" >> pipe_processing.log
echo "Specified species is $species" >> pipe_processing.log
echo "types of reads is $types" >> pipe_processing.log
echo " " >> pipe_processing.log

###################################################################################################
# Step 1, Trim ATAC-seq adapters and QC on seq file
# 1.1 Trim by cutadapt
if [[ $types == PE ]];
	then
	echo 'trimming ATAC PE reads by cutadapt'
	cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=15,10 --minimum-length=36  -o 'Trimed_'$name'_1.fastq' -p  'Trimed_'$name'_2.fastq'  $raw1 $raw2  > 'step1.1_'$name'_cutadapt_PE.trimlog'  
	temp=`grep "Total read pairs processed:" step1.1_*trimlog | awk '{print $5}'`  
	raw_reads=`echo ${temp//,}`  
elif [[ $types == SE ]];
	then
	echo 'trimming ATAC SE reads by cutadapt'
	cutadapt -a $adatper_1 --quality-cutoff=15,10 --minimum-length=36  -o  'Trimed_'$name'.fastq' $raw1 > 'step1.1_'$name'_cutadapt_SE.trimlog'  
	temp=`grep "Total reads processed:" step1.1_*trimlog | awk '{print $4}'`
	raw_reads=`echo ${temp//,}` 
fi

if [ $? == 0 ] 
	then
	echo "step1.1, trimming process sucessful!" >> pipe_processing.log
else 
	echo "step1.1, trimming process fail......" >> pipe_processing.log
	exit 1
fi



# 1.2 fastqc
echo 'fastqc is processing fastq file......'
fastqc -t $threads 'Trimed_'$name*'.fastq' -o . 

if [ $? == 0 ] 
	then
	echo "step1.2, fastqc process sucessful!" >> pipe_processing.log
else 
	echo "step1.2, fastqc process fail......" >> pipe_processing.log
	exit 1
fi

for zip in `ls | grep fastqc.zip`
do
unzip -o $zip	
mv $zip 'step1.2_'$zip
done


# 1.3 fastqc data collection (rely on the output data structure of Fastqc, double-check if it's updated)
echo -e "filename\tdeduplication_percentage\tmarker" > 'dedup_percentage_'$name'.result'
for file in `ls -d *fastqc/`
do
	cd $file
	temp=`echo ${file##Trimed_}`
	out_name=`echo ${temp%*_fastqc/}`
	out_value=`grep 'Total Deduplicated Percentage' fastqc_data.txt | awk '{print $4}'`
	echo -e "$out_name\t$out_value\t$marker" >> ../'dedup_percentage_'$name'.result'
	echo -e "item\t$out_name\t$out_name" > 'duplication_summary_'$out_name'.result'
	grep 'Sequence Duplication Levels' -A 15 fastqc_data.txt >> 'duplication_summary_'$out_name'.result'
	mv 'duplication_summary_'$out_name'.result' ../'data_collection_'$name
	echo -e "$out_name\tfastqc_test" > 'fastqc_summary_'$out_name'.result'
	awk -F "\t" '{print $1,$2}' OFS='\t' summary.txt >> 'fastqc_summary_'$out_name'.result'
	mv 'fastqc_summary_'$out_name'.result' ../'data_collection_'$name
	cd ..
done

if [ $? == 0 ] 
	then
	echo "step1.3, fastqc data_collection process sucessful!" >> pipe_processing.log
else 
	echo "step1.3, fastqc data_collection process fail......" >> pipe_processing.log
fi

mv 'dedup_percentage_'$name'.result'  ./'data_collection_'$name

# 1.4, get PE data R1 R2 deduplication difference percentage 
if [[ $types == PE ]];
then
per1=`tail -n 2 ./'data_collection_'$name/'dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '1p'`
per2=`tail -n 2 ./'data_collection_'$name/'dedup_percentage_'$name'.result' | awk '{print $2}' | sed -n '2p'`
dif=`echo "scale=2; ($per1-$per2)*200/($per1+$per2)" | bc -l`
else
dif=0
fi

if [ $? == 0 ] 
	then
	echo "step1.4, calculate replicate difference process sucessful!" >> pipe_processing.log
else 
	echo "step1.4, calculate replicate difference process fail......" >> pipe_processing.log
fi



###################################################################################################

# step2, alignment and data processing
# 2.1 alignment by bwa mem
echo 'alignment by bwa......'
bwa mem -t $threads $bwa_ref 'Trimed_'$name*'.fastq' | samtools view -bS - | samtools sort - -O 'bam' -o 'Trimed_'$name'.bam' -T temp_aln

if [ $? == 0 ] 
	then
	echo "step2.1, bwa alignment process sucessful!" >> pipe_processing.log
else 
	echo "step2.1, bwa alignment process fail......" >> pipe_processing.log
	exit 1
fi

# clean folder
find . -maxdepth 1 -name "Trimed*" ! -name "*bam" ! -name "step*" | xargs rm -r

for file in `ls *fastq 2> /dev/null`
do
if [[ $file != $R1 ]] && [[ $file != $R2 ]]
then
rm $file 2> /dev/null
fi
done

# 2.2, count distri after removing mapQ=0 (bwa would assign those reads to dif chr, and results in unreliable number)
samtools view -h 'Trimed_'$name'.bam' > input.sam
awk '$5>0' input.sam | sed '/@/d' - | cat <(grep '@' input.sam) - > output.sam
rm input.sam

# all alignment with mapQ > 0, exact same with samtools idxstats
cat output.sam | awk '{print $3}' |  sort  -k1,1V |   uniq -c > count_no_mapq0.txt
awk '! /random/ && ! /Un/ && /chr/  ' count_no_mapq0.txt  | awk '{print $2, $1}'  OFS="\t"  | sort  -k1,1 -V -s > temp.txt
if [[ $types == SE ]]; then
mv temp.txt  'chrom_count_'$name'.txt'
elif [[ $types == PE ]]; then
awk '{$2=int($2*0.5); print}' OFS="\t" temp.txt > 'chrom_count_'$name'.txt'
fi

# only effect reads
methylQA density -S $chrom_size  output.sam
cat output.extended.bed | awk '{print $1}' | uniq -c >  count_unique.txt
awk '! /random/ && ! /Un/ && /chr/  ' count_unique.txt  | awk '{print $2, $1}'  OFS="\t"  | sort  -k1,1 -V -s > 'chrom_count_unique_'$name'.txt'

# mapping status
map_mapped=`grep 'mappable reads' output.report | awk '{print $4}'`
map_uniq=`grep '(mapQ >= 10)' output.report | awk '{print $8}'`
map_effect=`grep 'non-redundant'  output.report | awk '{print $6}'`
mapped_ratio=`echo "scale=2; $map_mapped/$raw_reads" | bc -l`
effect_ratio=`echo "scale=2; $map_effect/$raw_reads" | bc -l`
nodup_ratio=`echo "scale=2; $map_effect/$map_uniq" | bc -l`

# rm chrM and other 
awk '$3!="chrM"' output.sam | samtools view -bS - > 'Trimed_rm_mapq0_chrm_'$name'.bam'
rm output*
rm count*.txt
awk -F "\t"  '{print $2}' 'chrom_count_unique_'$name'.txt'  |  paste 'chrom_count_'$name'.txt'  - | awk -F "\t" -v marker=$marker '{print $0, marker}' OFS="\t"  > ./'data_collection_'$name/'chrom_count_'$name'.result'

if [ $? == 0 ] 
	then
	echo "step2.2, count reads distribution process sucessful!" >> pipe_processing.log
else 
	echo "step2.2, count reads distribution process fail......" >> pipe_processing.log
fi

rm chrom_count*txt
mv 'Trimed_'$name'.bam'  'step2.1_Trimed_'$name'.bam'


###################################################################################################
# step3,QC and peak calling
# 3.1 QC by methylQA
echo 'methylQA processing......'
methylQA atac $chrom_size  'Trimed_rm_mapq0_chrm_'$name'.bam'

if [ $? == 0 ] 
	then
	echo "step3.1, mathylQA atac process sucessful!" >> pipe_processing.log
else 
	echo "step3.1, mathylQA atac process fail......" >> pipe_processing.log
	exit 1
fi

useful_reads=`grep 'non-redundant'  Trimed_rm_mapq0_chrm_*report | awk '{print $6}'`
mv 'Trimed_rm_mapq0_chrm_'$name'.bam'   'step2.2_Trimed_rm_mapq0_chrm_'$name'.bam'
mv 'Trimed_rm_mapq0_chrm_'$name'.genomeCov.pdf'  'step2.2_Trimed_rm_mapq0_chrm_'$name'.genomeCov'
awk '$1<=500'  'Trimed_'*$name'.insertdistro'  | sort -n | uniq -c | awk '{print $2,$1}' > 'insertion_distri_'$name'.result'
mv 'insertion_distri_'$name'.result'  ./'data_collection_'$name
rm 'Trimed_rm_mapq0_chrm_'$name'.insertdistro'*


# 3.2 normalization of *.bedGraph file by 10 Million reads
echo 'normalization bedGraph......'
norm=`grep 'non-redundant'  Trimed*report | awk '{print $6}'`
factor=`echo "scale=2; $norm/10000000" | bc -l`
awk -v factor=$factor '{print $1,$2,$3,$4/factor}' OFS='\t' 'Trimed_'*$name'.open.bedGraph'  >  'Normalized_per_10M_'$name'.open.bedGraph'
bedGraphToBigWig  'Normalized_per_10M_'$name'.open.bedGraph'   $chrom_size  'Normalized_per_10M_'$name'.bigWig'

if [ $? == 0 ] 
	then
	echo "step3.2, normalization process sucessful!" >> pipe_processing.log
else 
	echo "step3.2, normalization process fail......" >> pipe_processing.log
fi

mv 'Normalized_per_10M_'$name'.open.bedGraph'  'step3.2_Normalized_per_10M_'$name'.open.bedGraph'
mv 'Normalized_per_10M_'$name'.bigWig'  'step3.2_Normalized_per_10M_'$name'.bigWig'

# add a new bigwig file without black list
intersectBed -a 'Trimed_'*$name'.open.bedGraph'  -b $black_list -v > rmbl.bedGraph
bedGraphToBigWig rmbl.bedGraph $chrom_size 'step3.2_Trimed_nochrm_rmbl_'$name'.bigWig'
rm rmbl.bedGraph



# 3.3 peak calling
echo 'peak calling......'
awk '$2 < $3' 'Trimed_rm_mapq0_chrm_'$name'.open.bed' | awk  '{ if ((length($1) < 6) && (length($1) > 1))  print $0}' OFS='\t' > temp.open.bed
intersectBed  -a temp.open.bed  -b $black_list   -v > 'Trimed_rmbl_'$name'.open.bed'
rm temp.open.bed
rm 'Trimed_rm_mapq0_chrm_'$name'.open.bed'
macs2 callpeak -t ./'Trimed_rmbl_'$name'.open.bed' -g $macs2_genome -q 0.01 -n 'peakcall_'$name  --keep-dup 1000 --nomodel --shift 0 --extsize 150

if [ $? == 0 ] 
	then
	echo "step3.3, macs2 peak calling process sucessful!" >> pipe_processing.log
else 
	echo "step3.3, macs2 peak calling process fail......" >> pipe_processing.log
	exit 1
fi

# peak length distribution:
awk '{print $3-$2+1}' 'peakcall_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' > 'peak_length_distri_'$name'.result'
mv 'peak_length_distri_'$name'.result'  ./'data_collection_'$name

###################################################################################################
# step4, additional analysis
# 4.1 get reads under peak data
total=`wc -l 'Trimed_rmbl_'$name'.open.bed'|awk '{print $1}'`
echo "calculating reads under peak ratio......"
python $pipe_path'/try_compare.py'  data/   'Trimed_rmbl_'$name'.open.bed'    'peakcall_'$name'_peaks.narrowPeak'

if [ $? == 0 ] 
	then
	echo "step4.1, calculate reads under peak process sucessful!" >> pipe_processing.log
else 
	echo "step4.1, calculate reads under peak process fail......" >> pipe_processing.log
fi

sum=`awk '{s+=$5}END{print s}' reads.txt`
ratio=`echo "scale=2; $sum*100/$total" | bc -l`


# 4.2 enrichment ratio
# 4.2.1, top 20k peaks enrichment (it looks that this depends on depth)
noise_read=`intersectBed -a 'Trimed_rmbl_'$name'.open.bed'  -b 'peakcall_'$name'_peaks.narrowPeak' -f 0.5 -v | wc -l`  
denominator=`echo "scale=10; $noise_read / $genome_size " | bc -l`

sort -k1,1V -k2,2n reads.txt > sorted_read.txt
sort -k1,1V -k2,2n 'peakcall_'$name'_peaks.narrowPeak' | awk '{print $9}' | paste sorted_read.txt - | awk '{print $0}' OFS='\t' > rpkm_for_all_peak.Peak
sort  -k 7 -nr rpkm_for_all_peak.Peak  > sorted_rpkm_all_peak.Peak

peak_number=`wc -l sorted_rpkm_all_peak.Peak | awk '{print $1}'`

cal_enrich() {
rm top_peak.txt  2> /dev/null
rm top_ratio.txt  2> /dev/null
rm 'enrichment_'$name'.result'  2> /dev/null
start=$1
max=$2
for i in `seq $start 1000 $max`
do
	head -$i sorted_rpkm_all_peak.Peak > 'top_'$i'_peaks_by_qvalue.Peak'
	peak_length=`awk -F " " '{s+=($3-$2+1)}END{print s}' 'top_'$i'_peaks_by_qvalue.Peak'`
	read_number=`awk -F " " '{s+=$5}END{print s}' 'top_'$i'_peaks_by_qvalue.Peak'`
	echo $i >> top_peak.txt
	echo "scale=2; $read_number / $peak_length / $denominator" | bc -l >> top_ratio.txt
	rm 'top_'$i'_peaks_by_qvalue.Peak'
done
cut -f1 top_peak.txt  | paste -s | awk '{print "file",$0}' OFS="\t"  > enrich_peak.txt
cut -f1 top_ratio.txt  | paste -s | awk -v file=$name '{print file,$0}' OFS="\t"  > enrich_ratio.txt
cat enrich_peak.txt enrich_ratio.txt  > 'enrichment_'$name'.result'
}

if [ $peak_number -lt 1000 ]; then
	cal_enrich $peak_number $peak_number
elif [ $peak_number -ge 1000 ]; then
	bins=`echo "$peak_number / 1000" | bc`
	if [ $bins -ge 20 ]; then
	bins=20
	fi
	max=`echo "$bins * 1000" | bc`
	cal_enrich 1000 $max
fi

if [ $? == 0 ] 
	then
	echo "step4.2.1, calculate top20k enrichment ratio process sucessful!" >> pipe_processing.log
else 
	echo "step4.2.1, calculate top20k enrichment ratio process fail......" >> pipe_processing.log
fi

mv 'enrichment_'$name'.result'  ./'data_collection_'$name
rm *rpkm*
rm top*
rm enrich_*.txt

# 4.2.3, all peak enrichment (similar to 4.2.1, depends on depth)
# 4.2.2, coding promoter enrichment (better than previous 2)
# coding enrichment = ( reads in promoter / promoter length)  /  (total reads / genome size)
peak='peakcall_'$name'_peaks.narrowPeak'
bed='Trimed_rmbl_'$name'.open.bed'
echo "the bed file is $bed......"
echo "the peak file is $peak......"
total_reads=`wc -l $bed | awk '{print $1}'`
denominator=`echo "scale=10; $total_reads / $genome_size" | bc -l`
intersectBed -a $peak -b $coding_promoter -u > promoter_peak.bed
reads_in_promoter=`intersectBed -a $bed -b promoter_peak.bed -f 0.5 -u | wc -l | awk '{print $1}'`
promoter_number=`intersectBed -a $coding_promoter -b promoter_peak.bed -F 0.5 -u | wc -l | awk '{print $1}'`
promoter_length=`echo "$promoter_number * 2000" | bc -l`
enrichment_ratio=`echo "scale=3; $reads_in_promoter / $promoter_length / $denominator" | bc -l`
if [ $? == 0 ] 
	then
	echo "step4.2.2, coding promoter enrichment ratio process sucessful!" >> pipe_processing.log
else 
	echo "step4.2.2, coding promoter enrichment ratio process fail......" >> pipe_processing.log
fi
echo -e "name\ttotal_reads\tpromoter_number\treads_in_promoter\tenrichment_ratio" > 'enrichment_ratio_in_promoter_'$name'.result'
echo -e "$name\t$total_reads\t$promoter_number\t$reads_in_promoter\t$enrichment_ratio" >> 'enrichment_ratio_in_promoter_'$name'.result'
echo "the coding promoter enrichment for $name is $enrichment_ratio"
mv 'enrichment_ratio_in_promoter_'$name'.result'  'data_collection_'$name

# all peak enrichment
# (reads_under_peak/peak_length)/(read_outof_peak/(genome-peak_length))
reads_under_peak=`intersectBed -a $bed -b $peak -f 0.5 -u | wc -l`
peak_length=`awk '{s+=$3-$2+1}END{print s}' $peak`

all_peak_enrichment=`echo "scale=3; $reads_under_peak * ($genome_size - $peak_length)  / $peak_length / ($total_reads - $reads_under_peak)" | bc -l`
if [ $? == 0 ] 
	then
	echo "step4.2.3, all peak enrichment ratio process sucessful!" >> pipe_processing.log
else 
	echo "step4.2.3, all peak enrichment ratio process fail......" >> pipe_processing.log
fi
echo -e "name\ttotal_reads\tcoverage\treads_under_peak\tall_peak_enrichment" > 'all_peak_enrichment_'$name'.result'
echo -e "$name\t$total_reads\t$peak_length\t$reads_under_peak\t$all_peak_enrichment" >> 'all_peak_enrichment_'$name'.result'
echo "the enrichment for $name is $all_peak_enrichment"
mv 'all_peak_enrichment_'$name'.result' ./data_collection_*
unset peak bed


# 4.3 PBC calculation
methylQA atac -r -o PBC_calculation  $chrom_size  'step2.2_Trimed_rm_mapq0_chrm_'$name'.bam'
M_distinct=`awk '{s+=($3-$2+1)}END{print s}'  PBC_calculation.open.bedGraph`
M1=`awk '{if($4==1) s+=($3-$2+1)}END{print s}' PBC_calculation.open.bedGraph`
PBC=`echo "scale=2; $M1/$M_distinct" | bc -l`
rm PBC*
echo -e "file\ttotal\tmapped\tmapped_ratio\tuniq_mapped\tnon-redundant_uniq_mapped\teffect_ratio\tuseful_reads\tPBC\tnodup_ratio\tnumber_of_reads_under_peak\trup_ratio\treplicate_dif\tmarker" >  'mapping_status_'$name'.result'
echo -e "$name\t$raw_reads\t$map_mapped\t$mapped_ratio\t$map_uniq\t$map_effect\t$effect_ratio\t$useful_reads\t$PBC\t$nodup_ratio\t$sum\t$ratio\t$dif\t$marker" >>  'mapping_status_'$name'.result'
mv 'mapping_status_'$name'.result'  ./'data_collection_'$name

if [ -z $name ] || [ -z $raw_reads ] || [ -z $map_mapped ] || [ -z $mapped_ratio ] || [ -z $map_uniq ] || [ -z $useful_reads ] || [ -z $map_effect ] || [ -z $effect_ratio ] || [ -z $PBC ] || [ -z $nodup_ratio ] || [ -z $sum ]|| [ -z $ratio ]|| [ -z $dif ]
then
	echo "step4.3, sumarizing result process fail......" >> pipe_processing.log
else
	echo "step4.3, sumarizing result process sucessful!" >> pipe_processing.log
fi


# 4.4 IDR: Irreproducible Discovery Rate (IDR) based on statistical model
## 1), get 2 pseudoreplicates by random sampling bed file without replacement (original bed -> 2 new bed)
shuf 'Trimed_rmbl_'$name'.open.bed'  -o  random_reads.bed
n=`wc -l random_reads.bed | awk '{print $1}' `
size=$(( n/2 ))
head -$size random_reads.bed > sample1.bed
tail -$size random_reads.bed > sample2.bed
macs2 callpeak -t sample1.bed  -g $macs2_genome -q 0.01 -n   sample1_macs2 --keep-dup 1000 --nomodel --shift 0 --extsize 150
macs2 callpeak -t sample2.bed  -g $macs2_genome -q 0.01 -n   sample2_macs2 --keep-dup 1000 --nomodel --shift 0 --extsize 150

cp -r $idr_file  ./idr
mv sample*Peak  ./idr
cd ./idr
rm genome_table.txt
cp $chrom_size  ./genome_table.txt
Rscript batch-consistency-analysis.r  sample1_macs2_peaks.narrowPeak  sample2_macs2_peaks.narrowPeak  -1  $name'_self_IDR' 0 F p.value
if [ $? == 0 ] 
	then
	echo "step4.4, IDR process sucessful!" >> ../pipe_processing.log
else 
	echo "step4.4, IDR process fail......" >> ../pipe_processing.log
fi
Rscript batch-consistency-plot.r  1  plot_IDR  $name'_self_IDR' 
mkdir $name'_self_IDR'
mv $name'_self'*  ./$name'_self_IDR'  2> /dev/null
mv plot_IDR*  ./$name'_self_IDR'
mv $name'_self_IDR' ../
cd ..
rm -r ./idr
rm random_reads.bed
rm sample*
rm *zip
cp ./$name'_self_IDR'/plot_IDR-plot.ps  ./'data_collection_'$name/'idr_plot_'$name'.ps'


# 4.5 saturation analysis
# subsampling:
for number in 5 10 20 30 40 50 60 70 80 90
do
	sample_ratio=$(( total * $number / 100 ))
	shuf 'Trimed_rmbl_'$name'.open.bed' | head -$sample_ratio >  'Trimed_rmbl_'$name'_sample'$number'.open.bed'
done

if [ $? == 0 ] 
	then
	echo "step4.5, saturation subsampling process sucessful!" >> pipe_processing.log
else 
	echo "step4.5, saturation subsampling process fail......" >> pipe_processing.log
fi

# call peak
mkdir saturation_$name
mv *.open.bed ./saturation_$name/
cp peakcall_*Peak ./saturation_$name/'peakcall_Trimed_rmbl_'$name'.open.bed_peaks.narrowPeak'
cd ./saturation_$name
for file in `ls *sample*.open.bed`;
do
	macs2 callpeak -t $file -g $macs2_genome -q 0.01 -n 'peakcall_'$file --keep-dup 1000 --nomodel --shift 0 --extsize 150
done

if [ $? == 0 ] 
	then
	echo "step4.5, saturation call peak process sucessful!" >> ../pipe_processing.log
else 
	echo "step4.5, saturation call peak process fail......" >> ../pipe_processing.log
fi

echo "peak calling done......"

# summarise results
echo -e "5\n`seq 10 10 100`" > saturation_points.txt

for file in `ls *open.bed`
do
read_num=`wc -l $file | awk '{print $1}'`
echo `echo "scale=2; $read_num / 1000000" | bc -l`>> temp1.txt 
done
sort -k1,1n temp1.txt > saturation_reads.txt
rm temp1.txt

total_region=`awk '{s+=$3-$2+1}END{print s}' 'peakcall_Trimed_rmbl_'$name'.open.bed_peaks.narrowPeak'`

for file in `ls *narrowPeak`
do
peak_number=`wc -l $file | awk '{print $1}'`
peak_region=`intersectBed -a $file -b 'peakcall_Trimed_rmbl_'$name'.open.bed_peaks.narrowPeak' | awk '{s+=$3-$2+1}END{print s}'`
if [ -z "$peak_region" ]; then
peak_region=0
fi
echo `echo "scale=2; $peak_region / $total_region" | bc -l` >> temp3.txt
echo $peak_number >> temp2.txt
done

if [ $? == 0 ] 
	then
	echo "step4.5, saturation results collection process sucessful!" >> ../pipe_processing.log
else 
	echo "step4.5, saturation results collection process fail......" >> ../pipe_processing.log
fi


sort -k1,1n temp2.txt > saturation_peak.txt
sort -k1,1n temp3.txt > saturation_ratio.txt
rm temp2.txt
rm temp3.txt
paste saturation_points.txt  saturation_reads.txt  saturation_peak.txt   saturation_ratio.txt  > temp4.txt

echo -e "file\t$name'_read'\t$name'_peak'\t$name'_ratio'\tmarker"  > 'saturation_'$name'.result'  
awk -v marker=$marker '{print $0,marker}' OFS='\t' temp4.txt >> 'saturation_'$name'.result'  
rm temp4.txt
rm saturation*.txt
rm *sample*.open.bed
mv *.open.bed ../
mv 'saturation_'$name'.result'  ../'data_collection_'$name
cd ..

# 4.6 calculate background
peak='peakcall_'$name'_peaks.narrowPeak'
bed='Trimed_rmbl_'$name'.open.bed'

# signal part
intersectBed -a $peak -b $promoter_file -u | awk '{print $1"\t"$2"\t"$3"\t""1""\t"$9}' > promoter.narrowPeak
intersectBed -a $peak -b $promoter_file -v | awk '{print $1"\t"$2"\t"$3"\t""0""\t"$9}' > non-promoter.narrowPeak

echo -e "num_peaks_in_promoter\tnum_peaks_in_non-promoter\tnum_reads_in_promoter_peaks\tnum_reads_in_non-promoter_peaks" > 'promoter_percentage_'$name'.result'

peak1=`wc -l promoter.narrowPeak | awk '{print $1}'`
peak2=`wc -l non-promoter.narrowPeak | awk '{print $1}'`
read1=`intersectBed -a $bed -b promoter.narrowPeak -u -f 0.50 | wc -l`
read2=`intersectBed -a $bed -b non-promoter.narrowPeak -u -f 0.50 | wc -l`

echo -e "$peak1\t$peak2\t$read1\t$read2" >> 'promoter_percentage_'$name'.result'
sed -i 's/^-e //' 'promoter_percentage_'$name'.result'

cat promoter.narrowPeak non-promoter.narrowPeak | sort -k5 -n -r > top10k.narrowPeak
python $pipe_path'/promoter_bin.py' top10k.narrowPeak
rm promoter.narrowPeak
rm non-promoter.narrowPeak
rm top10k.narrowPeak

# background noise part
awk  '{ if ((length($1) < 6) && (length($1) > 1))  print $0}' OFS='\t' $chrom_size > temp.txt
python2.7  $pipe_path'/random_chr.py' temp.txt

size=`wc -l $bed | awk '{print $1}'`

awk '{print $1"\t"int(($3+$2)/2)-100000"\t"int(($3+$2)/2)+100000"\t"$4}' $peak > temp
awk '{if ($2<0) $2=0; print $0}' OFS="\t" temp > temp2
mv temp2 temp
intersectBed -a chr.peak -b temp -v | shuf - | head -50000 | sort -k1,1V -k2,2n > background
intersectBed -a $bed -b background -u -f 0.5 | sort -k1,1V -k2,2n > temp
python $pipe_path'/rpkm_bin.py' background temp $size

if [ $? == 0 ] 
	then
	echo "step4.6, background evaluation process sucessful!" >> ./pipe_processing.log
else 
	echo "step4.6, background evaluation process fail......" >> ./pipe_processing.log
fi

mv reads.txt 'background_'$name'.result'



rm temp
rm background
mv bin.txt 'bin_'$name'.result'

bg_total=`wc -l background*.result | awk '{print $1}'`  
bg_half_thres=`awk '$6<=0.15 {print $0}' background*.result | wc -l`  
bg_less=`awk '$6<=0.3 {print $0}' background*.result | wc -l`  
bg_more=`awk '$6>0.3 {print $0}' background*.result | wc -l` 
ra_half_thres=`echo "scale=2; $bg_half_thres*100 / $bg_total" | bc -l`  
ra_less=`echo "scale=2; $bg_less*100 / $bg_total" | bc -l`  
ra_more=`echo "scale=2; $bg_more*100 / $bg_total" | bc -l`  
echo -e "$ra_half_thres\t$ra_less\t$ra_more" > 'dichoto_bg_'$name'.result'  

mv 'dichoto_bg_'$name'.result'   ./'data_collection_'$name  
mv background*.result  ./'data_collection_'$name
mv promoter*.result ./'data_collection_'$name
mv bin*.result ./'data_collection_'$name

# step 4.7, plot on results
# clean result
find . -name "*.result" | xargs sed -i 's/^-e //'
cd ./'data_collection_'$name
Rscript $pipe_path'/visualization.R' $name
if [ $? == 0 ] 
	then
	echo "step4.7, plot process sucessful!" >> ../pipe_processing.log
else 
	echo "step4.7, plot process fail......" >> ../pipe_processing.log
fi

mkdir 'plots_collection_'$name
mv *png 'plots_collection_'$name
cp 'dedup_percentage_'$name'.result'   ../'step1.3_dedup_percentage_'$name'.result'
cp 'chrom_count_'$name'.result'  ../'step2.2_chrom_count_'$name'.result'
cp 'insertion_distri_'$name'.result'   ../'step3.1_insertion_distri_'$name'.result'
mv 'plots_collection_'$name  ../
cd ..

rm config.txt
rm -r data/
rm sorted_read.txt
rm chr.peak
rm -r 'saturation_'$name
rm temp.txt
rm refined_chrom_size.txt
rm pesudo_bl.txt 2> /dev/null

rename 's/Trimed_/step3.1_Trimed_/' Trimed_*
rename 's/peakcall_/step3.3_peakcall_/' peakcall_*

echo "Processing $name done"
echo "Processing $name done"
echo "Processing $name done"
cd ..
date







