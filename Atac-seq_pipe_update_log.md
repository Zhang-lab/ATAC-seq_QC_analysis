### Current version: v1.1a  
Last update: 02/25/2018  
  
  
  
### To be added soon:  
1, add pipe version and run time in output json file  
2, provide all results in ".txt" data format instead of figure directly in json output  
  
  
  
### Update record: 
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
 
