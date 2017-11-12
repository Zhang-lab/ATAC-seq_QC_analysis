import sys

chr=[]
with open(sys.argv[1],'r') as f:
	for line in f:
		line=line.strip('\n').split('\t')
		for i in range(int(line[1])/500):
			chr.append(line[0]+'\t'+str(i*500)+'\t'+str((i+1)*500)+'\t'+'random_non_peak_'+str(i+1)+'\n')
f.close()

with open('chr.peak','w') as f:
	f.writelines(chr)
f.close()
