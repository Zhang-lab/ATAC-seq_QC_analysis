### Current version: v1.2  
Last update: 04/11/2018  
  
  
  
### To be added soon:  

  
### Update record: 
04/11/2018, v1.2  
1. add warning for peak files with less than 100 peaks, there would be NO result for promoter percentage on peaks for those data  
2. reorganize pipe code for modulation  


03/26/2018, v1.1b  
1. add docker image id and md5sum for script file verification  
2. use soft link instead of mv for raw files (for cases in HTCF)  

02/25/2018, v1.1a  
Change output json file content. Put all raw data inside.  

02/23/2018  
Add insertion site record file in bedGraph and bigWig format for the purpose of narrowing down motif finding region  

02/06/2018  
Miner change in "chrom_count" result, add 0 for each cell so that there won't be error when some chrom has no count (e.g chrY=0).  

01/22/2018:  
1, Finalized version for pipe v1. Updated on server, Github and Docker (to be finished on the night of 01/22)  
2, Use subsample 10M calculation for enrichment instead of previous 40M normalization part, and update the reference for ENCODE PE data:  
enrichment=  
[ ($rupn+10000000*$peak_length / $genome_size) / $peak_length ] divided by	[ (10000000+$useful_ends) / ($genome_size-$peak_length)) ]  
3, change some hard coded variable names  
 
