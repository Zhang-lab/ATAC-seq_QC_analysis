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
# readlink or realpath can do the trick, but they are not available in some version, e.g. the -f option
pipe_loc=`dirname $0`
cd $pipe_loc
pipe_path=`pwd`

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
chrom_size="./refined_chrom_size.txt"

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
	cutadapt -a $adapter_1 -A $adapter_2 --quality-cutoff=10 --minimum-length=36  -o 'Trimed_'$name'_1.fastq' -p  'Trimed_'$name'_2.fastq'  $raw1 $raw2  > 'step1.1_'$name'_cutadapt_PE.trimlog'
elif [[ $types == SE ]];
	then
	echo 'trimming ATAC SE reads by cutadapt'
	cutadapt -a $adatper_1 --quality-cutoff=10 --minimum-length=36  -o  'Trimed_'$name'.fastq' $raw1 > 'step1.1_'$name'_cutadapt_SE.trimlog'
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
awk '! /random/ && ! /Un/ && /chr/  ' count_no_mapq0.txt  | awk '{print $2, $1}'  OFS="\t"  | sort  -k1,1 -V -s > 'chrom_count_'$name'.txt'

# only effect reads
methylQA density -S $chrom_size  output.sam
cat output.extended.bed | awk '{print $1}' | uniq -c >  count_unique.txt
awk '! /random/ && ! /Un/ && /chr/  ' count_unique.txt  | awk '{print $2, $1}'  OFS="\t"  | sort  -k1,1 -V -s > 'chrom_count_unique_'$name'.txt'

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
awk '{print $3-$2}' 'peakcall_'$name'_peaks.narrowPeak' | sort -n | uniq -c | awk '{print $2,$1}' > 'peak_length_distri_'$name'.result'
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


# add RUP results into mapping data collection
# mapping status
map_total=`grep "total reads" Trimed_*report | awk '{print $4}'`
map_mapped=`grep 'mappable reads' Trimed_*report | awk '{print $4}'`
map_uniq=`grep '(mapQ >= 10)' Trimed_*report | awk '{print $8}'`
map_effect=`grep 'non-redundant'  Trimed_*report | awk '{print $6}'`
mapped_ratio=`echo "scale=2; $map_mapped/$map_total" | bc -l`
uniq_ratio=`echo "scale=2; $map_uniq/$map_total" | bc -l`
effect_ratio=`echo "scale=2; $map_effect/$map_total" | bc -l`
nodup_ratio=`echo "scale=2; $map_effect/$map_uniq" | bc -l`


# 4.2 enrichment ratio
denominator=`echo "scale=10; $total / $genome_size " | bc -l`

sort -k1,1V -k2,2n reads.txt > sorted_read.txt
sort -k1,1V -k2,2n 'peakcall_'$name'_peaks.narrowPeak' | awk '{print $9}' | paste sorted_read.txt - | awk '{print $0}' OFS='\t' > rpkm_for_all_peak.Peak
sort  -k 7 -nr rpkm_for_all_peak.Peak  > sorted_rpkm_all_peak.Peak

peak_number=`wc -l sorted_rpkm_all_peak.Peak` | awk '{print $1}'

