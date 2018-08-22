#!/bin/bash
date_now=`date +"%m-%d-%y"`
mkdir 'result_summary_'$date_now
find . -name "QC_ATAC_data_collection_*"  -type d | xargs cp -r -t 'result_summary_'$date_now
cd 'result_summary_'$date_now


# 1, enrichment
# 1.1 all peak enrichment
echo -e "name\trup_reads\trup\tcoverage\tsub10M_peak_enrichment" > merged_sub10M_enrichment.txt
for file in `find . -name "step4.2_sub10M_enrichment*result"`
do
    sed -n '2p' $file >> merged_sub10M_enrichment.txt
    rm $file
done

# 1.2, coding promoter peak enrichment
echo -e "name\ttotal_reads\tpromoter_number\treads_in_promoter_peak\tenrichment_ratio" > merged_coding_promoter_peak_enrichment.txt
for file in `find . -name "step4.2_enrichment_ratio_in_promoter*.result"`
do
    sed -n '2p' $file >> merged_coding_promoter_peak_enrichment.txt
    rm  $file
done


# 2, mapping status
find . -name "QC_data_collection_*result" | xargs mv -t .
sed -n '1p'  `ls QC_data_collection_*result | head -1` > merged_mapping_status.txt
for file in `ls QC_data_collection_*result`
do
sed -n '2p' $file >> merged_mapping_status.txt
rm $file
done

find . -name "QC_table_*.result" | xargs mv -t .  2> /dev/null

if ls QC_table_*.result 1> /dev/null 2>&1; then
    sed -n '1p'  `ls QC_table_*.result | head -1` > merged_QC_table.txt
    for file in `ls QC_table_*.result`
	do
	sed -n '2p' $file >> merged_QC_table.txt
	rm $file
    done
fi

[ -f merged_QC_table.txt ] && rm merged_mapping_status.txt

# 2.2, unique_chrM
find . -name "step2.2_unique_chrM_ratio*result" | xargs mv -t .
echo -e "name\tunique_mapped_reads\tunique_chrM\tratio" > merged_unique_chrM_ratio.txt
for file in `ls step2.2_unique_chrM_ratio*result`
do
    temp=`echo ${file##step2.2_unique_chrM_ratio_}`
    name=`echo ${temp%%.result}`
    echo -e "$name\t`sed -n '2p' $file`" >> merged_unique_chrM_ratio.txt
    rm $file
done

# 2.3, useful reads
find . -name "step3.1_useful_reads*result" | xargs mv -t .
sed -n '1p'  `ls step3.1_useful_reads*result | head -1` > merged_useful_reads.txt
for file in `ls step3.1_useful_reads*result`
do
sed -n '2p' $file >> merged_useful_reads.txt
rm $file
done




# 3, background
# 3.1, random bg
find . -name "step4.5_background_*.result" | xargs mv -t .
for file in `ls step4.5_background_*.result`
do
temp=`echo ${file##step4.5_background_}`
name=`echo ${temp%%.result}`
awk '{print $6}' $file > temp.txt
echo "$name" | cat - temp.txt  > $file
done
paste step4.5_background_*result > merged_background.txt
rm step4.5_background_*result temp.txt

# 3.2, dichoto count
find . -name "step4.5_dichoto_bg_*.result" | xargs mv -t .
echo -e "file\tlt_0.188\tlt_0.377\tgt_0.377" > merged_bg_dichoto.txt
for file in `ls step4.5_dichoto_bg_*.result`
do
    temp=`echo ${file#step4.5_dichoto_bg_}`
    name=`echo ${temp%.result}`
    echo $name
    add=`cat $file`
    echo -e "$name\t$add" >> merged_bg_dichoto.txt
    rm $file
done


# 4, promoter status
# 4.1, promoter percentage
find . -name "step4.5_promoter_percentage_*.result" | xargs mv -t  .
for file in `ls step4.5_promoter_percentage_*.result`
do
temp=`echo ${file##step4.5_promoter_percentage_}`
name=`echo ${temp%%.result}`
awk -v name=$name '{if(NR==1) print "file",$0; if(NR==2) print name,$0}' OFS='\t' $file > temp.txt
mv temp.txt $file
done
cat step4.5_promoter_percentage_*.result  | sort | uniq > merged_promoter_percentage.txt
rm step4.5_promoter_percentage_*.result

# 4.2 100bin count
find . -name "step4.5_bin_*.result" | xargs mv -t  .
for file in `ls step4.5_bin_*.result`
do
temp=`echo ${file##step4.5_bin_}`
name=`echo ${temp%%.result}`
echo -e "bin_seq\t$name" > temp.txt
sed -i 's/^-e //' temp.txt
cat temp.txt $file  > temp2.txt
mv temp2.txt $file
done
paste step4.5_bin_*result  > temp.txt
awk '{line=$1;for(i=2;i<=NF;i+=2) line=line"\t"$i; print line}' temp.txt >  merged_bin_result.txt
rm step4.5_bin_*.result
rm temp.txt

# 5, deduplication percentage
find . -name "step1.3_dedup_percentage*.result" | xargs mv -t .
sed -n '1p'  `ls step1.3_dedup_percentage*.result | head -1` > merged_dedup_percentage.txt
for file in `ls step1.3_dedup_percentage*.result`
do
sed -n '2p' $file >> merged_dedup_percentage.txt
rm $file
done

# 6, by python:
# 6.1, chrom count
find . -name "step2.2_chrom_count*result" | xargs mv -t .
#6.2, dedup summary
find . -name "step1.3_duplication_summary*result" | xargs mv -t .
#6.3, fastqc
find . -name "step1.3_fastqc_summary*result" | xargs mv -t .
#6.4, saturation
mkdir saturation
find . -name "step4.4_saturation*result" | xargs mv -t saturation/

cat <<EOF > temp.py
# python code to merge
import pandas as pd
import glob
import os
import sys

#1, merge all qc summary
qc_list=glob.glob("step1.3_fastqc_summary*result")
total=[]
for i in qc_list:
    content=pd.read_csv(i, header=0, sep='\t')
    total.append(content)

result=total[0]
for i in range(len(total)-1):
    result=pd.merge(result, total[i+1], left_on='fastqc_test', right_on='fastqc_test', how='outer')

result.to_csv("merged_qc_summary.txt", sep="\t", header=True, index=False)


#2, merge all chrom count and generate chrM ratio
chrom_list=glob.glob("step2.2_chrom_count*result")
total=[]
ratio=[]

for i in chrom_list:
    content=pd.read_csv(i, sep='\t', header=None)
    name=i[12:-7]
    content.columns=['chrom', '%s_mapped' %name, '%s_effect' %name, '%s_marker' %name ]
    mapped=content['%s_mapped' %name].sum()/100
    effected=content['%s_effect' %name].sum()/100
    content['%s_mapped_distri' %name]=content['%s_mapped' %name].divide(mapped, axis=0)
    content['%s_effect_distri' %name]=content['%s_effect' %name].divide(effected, axis=0)
    content=content[['chrom', '%s_marker' %name , '%s_mapped' %name, '%s_mapped_distri' %name, '%s_effect' %name, '%s_effect_distri' %name]]
    total.append(content)


result=total[0]
for i in range(len(total)-1):
    result=pd.merge(result, total[i+1], on=['chrom'])

# merged data
result.to_csv("merged_chrom_count.txt", sep="\t", header=True, index=False)


#3, merge duplication summary
dup_list=glob.glob("step1.3_duplication_summary*result")
total=[]
for i in dup_list:
    content=pd.read_csv(i, sep='\t', header=0)
    total.append(content)

result=total[0]

for i in range(len(total)-1):
    result=pd.merge(result, total[i+1], on=['item'])

result.to_csv("merged_dup_level.txt", sep="\t", header=True, index=False)


#4, saturation
os.chdir("./saturation")
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import matplotlib.pyplot as plt
file=glob.glob("step4.4_saturation*result")
peak=[]
read=[]
ratio=[]
label=[]
for i in file:
    table=pd.read_csv(i, sep='\t', header=0)
    table.columns=['sample', 'read', 'peak', 'ratio', 'marker']
    table=table.loc[:, ['read', 'peak','ratio']]
    table=pd.DataFrame(table, dtype='float')
    read.append(table['read'])
    peak.append(table['peak'])
    ratio.append(table['ratio'])
    name=i[11:-7]
    label.append(name)

for i in range(0, len(file)):
    plt.plot(read[i],peak[i], label=label[i])

plt.legend(loc='upper left', bbox_to_anchor=(-0.15, 1.1), ncol=1, fancybox=True, shadow=True, fontsize=10)
plt.xlabel("million_reads")
plt.ylabel("peak_numbers")
plt.savefig("peak_by_reads_mark.png")
plt.close()

for i in range(0, len(file)):
    plt.plot(read[i],ratio[i], label=label[i])

plt.legend(loc='upper left', bbox_to_anchor=(-0.15, 1.1), ncol=1, fancybox=True, shadow=True, fontsize=10)
plt.ylim=([0,1])
plt.axhline(y=0.8)
plt.xlabel("Million_reads")
plt.ylabel("ratio_of_peaks")
plt.savefig("ratio_by_reads.png")
plt.close()

# store merged file
table=pd.DataFrame()
for i in range(0, len(file)):
    table["%s_read" %label[i]]=read[i]
    table["%s_peak" %label[i]]=peak[i]
    table["%s_ratio" %label[i]]=ratio[i]

table.to_csv("merged_saturation_collection.txt", sep='\t', header=True)
EOF

python3.5  temp.py && rm temp.py
rm *result
cp ./saturation/merged_saturation_collection.txt  .

mkdir peak_length_distri
find . -name "step3.4_peak_length*result" | xargs mv -t ./peak_length_distri

mkdir insertion_size
find . -name "step3.1_insertion*result" | xargs mv -t ./insertion_size

cd saturation
mkdir single_record
mv step4.4_saturation*result  ./single_record
cd ..

mkdir preseq_estimate
find . -name "step2.3_yield_*result" | xargs mv -t  ./preseq_estimate

rm -r QC_ATAC_data_collection_*
cd ..

date
echo "process finished......."












