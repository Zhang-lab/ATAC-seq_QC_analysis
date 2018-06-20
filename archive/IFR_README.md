# Insertion Free Region Identification

06/18/2018
## Wellington Algorithm
Wellington algorithm is a novel method for the accurate identification of digital genomic footprints from DNase-seq data 
proposed by [*Jason Piper etc. 2013*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3834841/), which shares the same 
charateristics with identifiction of insertion free region within open chromatins from ATAC-seq data. The 
[`ATAC-seq_wellington.R`](../pipe_code/ATAC-seq_wellington.R) is an implemention of applying Wellington algorithm to peaks 
from ATAC-seq data called by MACS2.

Using the same notations of the original paper, let <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}"/> be 
the length (in base pairs) of the possible IFR, and <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/> be 
the length (in base pairs) of the shoulder on each side of the possible IFR. We consider three insertion counts in these 
regions: the total number of insetions inside the possible IFR 
(<img src="https://latex.codecogs.com/svg.latex?\Large&space;{FP}"/>), the insertion count in the upstream shoulder region on 
the forward reference strand (<img src="https://latex.codecogs.com/svg.latex?\Large&space;{SH}^+"/>), and the insertion count 
in the downstream shoulder region on the backward reference strand 
(<img src="https://latex.codecogs.com/svg.latex?\Large&space;{SH}^-"/>). 

We will test the null hypothesis that the number of insertions is proportional to the region length by using a binomial 
test. Considered the characteristics of our PEasSE peakcalling strategy, we will test the both strands together without losing 
sensitivity but having higher proficiency. With <img src="https://latex.codecogs.com/svg.latex?\Large&space;F(k,n,p)"/> being 
the binomial culmulative distribution function, i.e. the probability of achieving at least 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;k"/> out of 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;n"/> success with the probability of each success 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;p"/>, the *p*-value will be calculated by
<img src="https://latex.codecogs.com/svg.latex?\Large&space;F(FP,FP+{SH}^{+}+{SH}^-,l_{FP}/(l_{FP}+l_{SH}))"/>, which is given 
for a combination of possible <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}"/> and 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/>.

Then, *p*-values for all posisble combination of 
different <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}"/> and 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/> will be calculated independently and adjusted for 
multiple testing by the user-specified method. We can identify IFR within the peaks by using an appropriate threshold for the 
*p*-value and subsequently using a greedy selection strategy for IFR identification. The parameters 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}"/> and 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/> are individually determined for each IFR using maximum 
likelihood estimation. The combination of <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}"/> and 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/>, which leads to the smallest *p*-value as well as smaller 
than the given threshold, will be condsidered as the true possible IFR.

For a given peak, we first identified the most significant IFR within the whole peak. If there is an IFR, we then remove the 
identified IFR from the peak and search any other possible IFR within the rest part of the peak, until no IFR can be identified.

## Parameters Description
We have incorporated the IFR identication into our ATAC-seq pipeline, and user can simply add `-i` in the bash command to 
enable IFR identification with default parameters. User can also specify their own parameters by using 
`-i <min_lfp,max_lfp,step_lfp,min_lsh,max_lsh,step_lsh,method,p_cutoff>`. **These parameters are positional!!!**

(1) `min_lfp`: the minimum value of <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}/2"/> (default `5`).  
(2) `max_lfp`: the maximum value of <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}/2"/> (default `15`).  
(3) `step_lfp`: the increasing searching step of <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{FP}"/> from minimum 
to maximum (default `2`).  
(4) `min_lsh`: the minimum value of <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/> (default `50`).  
(5) `max_lsh`: the maximum value of <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/> (default `200`).  
(6) `step_lsh`: the increasing searching step of <img src="https://latex.codecogs.com/svg.latex?\Large&space;l_{SH}"/> from minimum 
to maximum (default `20`).  
(7) `method`: the method used for multiple testing adjustment (default `BH`). Possible values include `BH`, `BY`, `fdr`, `holm`, 
`bonferroni`, `hochberg` and `hommel`.  
(8) `p_cutoff`: the threshold of *p*-value (default `0.05`).

## Parallel Processing
Parallel processing is enabled by [GNU](https://www.gnu.org/software/parallel/parallel_tutorial.html), and shared the same 
number of thread in the whole pipeline process.
