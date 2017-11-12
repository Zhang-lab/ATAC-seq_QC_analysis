#!/bin/bash 
# This is for Docker pipe, default root dir is /atac-seq
# This is for Docker!!!
# This is for Docker!!!


species=$1

# get pipe path, though readlink/realpath can do it, some version doesn't have that
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
 

# common tools:
idr_file=$pipe_path'/../Resource/idr_folder'
adapter_1="CTGTCTCTTATACACATCT"
adapter_2="CTGTCTCTTATACACATCT"
fastq_dump_tool='fastq-dump.2.8.2'
preseq="preseq"

# genome specific resources:
if [[ $species == mm10 ]]; 
	then
	bwa_ref="/atac_seq/Resource/Genome/mm10/bwa_index_mm10/mm10.fa"
	chrom_size="/atac_seq/Resource/Genome/mm10/mm10.chrom.sizes"
	black_list="/atac_seq/Resource/Genome/mm10/mm10_black_list.bed"
	genome_size=2730871774
	promoter_file="/atac_seq/Resource/Genome/mm10/mm10_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/mm10/mm10_promoter_coding_bistream_1kb.bed"
	macs2_genome='mm'
elif [[ $species == mm9 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/mm9/bwa_index_mm9/mm9.fa"
	chrom_size="/atac_seq/Resource/Genome/mm9/mm9_chrom_sizes"
	black_list="/atac_seq/Resource/Genome/mm9/mm9-blacklist.bed"
	genome_size=2725765481
	promoter_file="/atac_seq/Resource/Genome/mm9/mm9_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/mm9/mm9_coding_promoter_bistream_1kb.bed"
	macs2_genome='mm'
elif [[ $species == hg38 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/hg38/bwa_index_hg38.25/hg38.25_chromsome.fa"
	chrom_size="/atac_seq/Resource/Genome/hg38/hg38.25_chromsome.sizes"
	black_list="/atac_seq/Resource/Genome/hg38/hg38_black_list.bed"
	genome_size=3209286105
	promoter_file="/atac_seq/Resource/Genome/hg38/hg38_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/hg38/hg38_coding_promoter_bistream_1kb.bed"
	macs2_genome='hs'
elif [[ $species == hg19 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/hg19/bwa_index_0.7.5/hg19.fa"
	chrom_size="/atac_seq/Resource/Genome/hg19/hg19_chromosome.size"
	black_list="/atac_seq/Resource/Genome/hg19/hg19_blacklist.bed"
	genome_size=3137161264
	promoter_file="/atac_seq/Resource/Genome/hg19/hg19_promoter_bistream_1kb.bed"
	coding_promoter="/atac_seq/Resource/Genome/hg19/hg19_promoter_coding_bistream_1kb.bed"
	macs2_genome='hs'
elif [[ $species == danRer10 ]];
	then
	bwa_ref="/atac_seq/Resource/Genome/danRer10/bwa_index_denRer10/danRer10.fa"
	chrom_size="/atac_seq/Resource/Genome/danRer10/danRer10.chrom.sizes"
	touch pesudo_bl.txt
	black_list="pesudo_bl.txt"
	genome_size=1371719383
	promoter_file="/atac_seq/Resource/Genome/danRer10/promoter_region_danRer10_bistream_1k.bed"
	coding_promoter="/atac_seq/Resource/Genome/danRer10/danRer10_coding_promoter_bistream_1k.bed"
	macs2_genome='mm'
elif [[ $species == personalize ]];
	then
	echo "please add all your preferred file as reference, please make sure that you are very clear of which file is for which"
	echo "remove exit 1 after adding your file"
	exit 1
	baw_ref=" "
	chrom_size=" "
	black_list=" "
	genome_size=
	promoter_file=" "
	macs2_genome=" "
	coding_promoter=" "
fi

# other tools note (their names are stable, so I use it directly)
# 1, cutadapt v1.12
# 2, fastqc v0.11.5
# 3, bwa v0.7.16a
# 4, methylQA v0.1.9
# 5, R version > v3.2.3
