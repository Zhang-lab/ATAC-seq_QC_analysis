
import pandas as pd
import sys
import numpy as np
import os

name=sys.argv[1]  # 'xxx_peaks.narrowPeak'
prefix=name[: -26]


point = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]
original=pd.read_csv('%s' %name, sep='\t', header=None)
total_peak = len(original)

ratio=[]
peak=[]
read=[]

for p in point:
   if os.stat("%s_sample%d.open.bed_peaks.narrowPeak" %(prefix, p)).st_size < 100:
      peaks=0
   else:
      table=pd.read_csv('%s_sample%d.open.bed_peaks.narrowPeak' %(prefix, p), sep='\t', header=None)
      peaks=len(table)
   peak.append(peaks)
   ratio.append(format(peaks / total_peak, '.2f'))

ratio.append(1)
point.append(100)
peak.append(total_peak)

prefix2=prefix[9:]
total_reads=len(pd.read_csv('%s.open.bed' %prefix2, sep='\t', header=None))
reads=round(total_reads/1000000, 2)
read[:]=[x*reads/100 for x in point]

data=[point,read,peak,ratio]
result=pd.DataFrame(data).transpose()
result.columns=['point','read','peak','ratio']
result.to_csv("%s_saturation_report.txt" %prefix2, header=True, sep='\t', index=False)

import matplotlib
matplotlib.use('Agg')

import pylab as pl
pl.plot(read, peak)
pl.plot(read, peak, 'ro')
pl.xlabel('reads_by_million')
pl.ylabel('peaks_number')
pl.axhline(y=peak[-1]*0.8)
pl.axhline(y=peak[-1])
pl.savefig("%s_saturation_plot_%d_peaks_total.png" %(prefix, total_peak))
