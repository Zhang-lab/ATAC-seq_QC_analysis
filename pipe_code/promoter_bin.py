import sys

peak=[]
with open(sys.argv[1],'r') as f:
	for line in f:
		line=line.strip('\n').split('\t')
		peak.append(int(line[3]))
f.close()

num=int(len(peak)/100.0)
bin=[]
for i in range(99):
	bin.append(str(i+1)+'\t'+str(sum(peak[num*i:num*(i+1)])/(num*1.0))+'\n')
bin.append('100'+'\t'+str(sum(peak[num*99:])/(num*1.0))+'\n')

with open('bin.txt','w') as f:
	f.writelines(bin)
f.close
