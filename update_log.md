### Current version: v4  
Last update: 08/23/2018  
  
  
### To be added soon:  

  
### Update record:  
`08/23/2018, v4`  
1. remove chrM reads in "non-redundant-uniquely-mapped reads"  
2. set default read length cutoff in methylQA to 38 instead of 50, use option -c to specify other numbers  
3. change unique chrM ratio formula, but resuit is same (now no chrM in report, line 336)  
4. change nodup_ratio formula, but result is same (now no chrM in report, line 343)  

`07/06/2018, v3.1`  
!!! Enrichment is overestimated due to chrM reads  
1. modify QC table: effect reads -> useful single ends; 2 dedup -> dup rate    
2. edit help information  
3. adjust saturation calculation: use overlapped region coverage  
4. correct visualization.R percentage calculation  
5. add R package "data.table" into docker  
6. fix that all "percentage" in json output are numeric values (0.01 for 1%)  

`06/28/2018, v3`  
1. add insertion free region finding algorithm  
2. remove parallel running  
3. modify QC table and output structure

`05/07/2018, targetv2` 
1. fix version for target (v2)  
2. fix docker image and image id for now  

`04/20/2018, v1.2b`  
1. modify `intersectBed` cmd for HTCF using, add `-iobuf 200M` for all of them   


`04/11/2018, v1.2`  
1. add warning for peak files with less than 100 peaks, there would be NO result for promoter percentage on peaks for those data  
2. reorganize pipe code for modulation  



`03/26/2018, v1.1b`  
1. add docker image id and md5sum for script file verification  
2. use soft link instead of mv for raw files (for cases in HTCF)  


`02/25/2018, v1.1a`  
Change output json file content. Put all raw data inside.  


`02/23/2018`  
Add insertion site record file in bedGraph and bigWig format for the purpose of narrowing down motif finding region  


`02/06/2018`  
Miner change in "chrom_count" result, add 0 for each cell so that there won't be error when some chrom has no count (e.g chrY=0).  


`01/22/2018`  
1, Finalized version for pipe v1. Updated on server, Github and Docker (to be finished on the night of 01/22)  
2, Use subsample 10M calculation for enrichment instead of previous 40M normalization part, and update the reference for ENCODE PE data:  
enrichment=  
[ ($rupn+10000000*$peak_length / $genome_size) / $peak_length ] divided by	[ (10000000+$useful_ends) / ($genome_size-$peak_length)) ]  
3, change some hard coded variable names  
 
