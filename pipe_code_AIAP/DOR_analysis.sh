#!/bin/bash
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#group should be the character that distinguish 2 bed/peak file
group1=$1
group2=$2

if [ -z "$group2" ]
    then
    echo "missing 2 group input!!!"
    echo "group should be the character that distinguish 2 bed/peak file"
    exit
fi

cat *peaks.narrowPeak | sort -k1,1V -k2,2n | bedtools merge -i - -c 5 -o max | awk '{print $1"\t"$2"\t"$3"\t""merged_peak"NR"\t"$4"\t""."}' > merged_all.narrowPeak
cut -f1-3 merged_all.narrowPeak > full.count

for open in `ls step3.3*_rmbl*open.bed | grep $group1`
do
    intersectBed -a $open -b merged_all.narrowPeak -u -f 0.5 | sort -k1,1V -k2,2n > temp
	python ${pipe_path}/rpkm_bin.py merged_all.narrowPeak temp 10
	paste full.count <(awk '{print $5}' reads.txt) > temp
	mv temp full.count
done

for open in `ls step3.3*_rmbl*open.bed | grep $group2`
do
    intersectBed -a $open -b merged_all.narrowPeak -u -f 0.5 | sort -k1,1V -k2,2n > temp
	python ${pipe_path}/rpkm_bin.py merged_all.narrowPeak temp 10
	paste full.count <(awk '{print $5}' reads.txt) > temp
	mv temp full.count
done

num1=`ls step3.3*rmbl*open.bed | grep $group1 | wc -l`
num2=`ls step3.3*rmbl*open.bed | grep $group2 | wc -l`

Rscript ${pipe_path}/DOR_analysis.R $group1 $group2 $num1 $num2
