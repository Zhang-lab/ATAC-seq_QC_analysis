#!/bin/bash  
species=$1

# get pipe path, though readlink/realpath can do it, some version doesn't have that
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
 

# executables and variables
adapter_1="CTGTCTCTTATACACATCT"
adapter_2="CTGTCTCTTATACACATCT"
cutadapt="/usr/bin/cutadapt"
fastqc="/usr/bin/fastqc"
preseq="/usr/bin/preseq"
fastq_dump_tool="/usr/bin/fastq-dump.2.8.1-3"
samtools="/usr/local/bin/samtools"
bwa="/usr/bin/bwa"
methylQA="/usr/bin/methylQA"
macs2="/usr/local/bin/macs2"

# genome specific resources:
if [[ $species == mm10 ]]; 
    then
    bwa_ref="/home/Resource/Genome/mm10/bwa_index_mm10/mm10.fa"
    chrom_size="/home/Resource/Genome/mm10/mm10.chrom.sizes"
    black_list="/home/shaopengliu/resources/black_list/mm10_black_list.bed"
    genome_size=2730871774
    promoter_file="/home/shaopengliu/resources/promoter_region/mm10_promoter_bistream_2kb.bed"
    coding_promoter="/home/shaopengliu/resources/promoter_region/mm10_promoter_coding_bistream_1kb.bed"
    macs2_genome='mm'
elif [[ $species == mm9 ]];
    then
    bwa_ref="/home/Resource/Genome/mm9/bwa_index_mm9/mm9.fa"
    chrom_size="/home/Resource/Genome/mm9/mm9_chrom_sizes"
    black_list="/home/Resource/Genome/mm9/mm9-blacklist.bed"
    genome_size=2725765481
    promoter_file="/home/shaopengliu/resources/mm9/mm9_promoter_bystream_2kb.bed"
    coding_promoter="/home/shaopengliu/resources/mm9/coding_promoter_mm9.bed"
    macs2_genome='mm'
elif [[ $species == hg38 ]];
    then
    bwa_ref="/home/Resource/Genome/hg38/bwa_index_hg38.25/hg38.25_chromsome.fa"
    chrom_size="/home/Resource/Genome/hg38/hg38.25_chromsome.sizes"
    black_list="/home/shaopengliu/resources/black_list/hg38_black_list.bed"
    genome_size=3209286105
    promoter_file="/home/shaopengliu/resources/promoter_region/hg38_promoter_bistream_2kb.bed"
    coding_promoter="/home/shaopengliu/resources/promoter_region/hg38_coding_promoter_bistream_1kb.bed"
    macs2_genome='hs'
elif [[ $species == hg19 ]];
    then
    bwa_ref="/home/Resource/Genome/hg19/bwa_index_0.7.5/hg19.fa"
    chrom_size="/home/Resource/Genome/hg19/hg19_chromosome.size"
    black_list="/home/Resource/Genome/hg19/hg19_blacklist.bed"
    genome_size=3137161264
    promoter_file="/home/shaopengliu/resources/hg19/hg19_promoter_bistream_1kb.bed"
    coding_promoter="/home/shaopengliu/resources/hg19/hg19_promoter_coding_bystream_1kb.bed"
    macs2_genome='hs'
elif [[ $species == danRer10 ]];
    then
    bwa_ref="/home/Resource/Genome/danRer10/bwa_index_denRer10/danRer10.fa"
    chrom_size="/home/shaopengliu/resources/danRer10/danRer10.chrom.sizes"
    touch pesudo_bl.txt
    black_list="pesudo_bl.txt"
    genome_size=1340447187
    promoter_file="/home/shaopengliu/resources/danRer10/promoter_region_danRer10_bistream_2k.bed"
    coding_promoter="/home/shaopengliu/resources/danRer10/danRer10_coding_promoter_bistream_2k.bed"
    macs2_genome='mm'
elif [[ $species == dm6 ]];
    then
    bwa_ref="/home/shaopengliu/resources/dm6/bwa_index_dm6/dm6.fa"
    chrom_size="/home/shaopengliu/resources/dm6/d.mel.chrom.sizes"
    touch pesudo_bl.txt
    black_list="pesudo_bl.txt"
    genome_size=143726002
    promoter_file="/home/shaopengliu/resources/dm6/promoter_region_from_Dmel.bed"
    coding_promoter="/home/shaopengliu/resources/dm6/pseudo_coding_promoter_region.bed"
    macs2_genome='dm'
elif [[ $species == rn6 ]];
    then
    bwa_ref="/home/Resource/Genome/rn6/bwa_index_rn6/rn6.fa"
    chrom_size="/home/Resource/Genome/rn6/rn6.chrom.sizes"
    touch pesudo_bl.txt
    black_list="pesudo_bl.txt"
    genome_size=2870182909
    promoter_file="/home/Resource/Genome/rn6/promoter_region.bed"
    coding_promoter="/home/Resource/Genome/rn6/coding_promoter_region.bed"
    macs2_genome='mm'
elif [[ $species == danRer11 ]];
    then
    bwa_ref="/home/Resource/Genome/GRCz11/bwa_index_GTCz11/GRCz11.fa"
    chrom_size="/home/Resource/Genome/GRCz11/GRCz11_chrom.size"
    touch pesudo_bl.txt
    black_list="pesudo_bl.txt"
    genome_size=1345118429
    promoter_file=/home/shaopengliu/resources/danRer11/GRCz11_promoter_region.bed
    coding_promoter=/home/shaopengliu/resources/danRer11/pseudo_GRCz11_promoter_region.bed    macs2_genome='mm'
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


