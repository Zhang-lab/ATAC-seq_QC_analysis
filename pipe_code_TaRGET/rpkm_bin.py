import sys

peak=[]
with open(sys.argv[1],'r') as f:
	for line in f:
		line=line.strip('\n').split('\t')
		peak.append(line)

bed=[]
with open(sys.argv[2],'r') as f:
	for line in f:
		line=line.strip('\n').split('\t')
		bed.append(line)

SIZE=int(sys.argv[3])

index=0
n=len(peak)
num=[0]*n
for read in bed:
	mid=(int(read[1])+int(read[2]))/2
	while (index<n-1 and mid>int(peak[index][2])) or (index<n-1 and read[0]!=peak[index][0]):
		index+=1
	num[index]+=1
	if (index<n-1) and (mid==int(peak[index+1][1])):
		num[index+1]+=1

output=[]
for i in range(n):
	if num[i]!=0:
		y=1.0*num[i]*10**9/SIZE/(int(peak[i][2])-int(peak[i][1]))
		y='%.4f'%y
		output.append(peak[i][0]+'\t'+peak[i][1]+'\t'+peak[i][2]+'\t'+peak[i][3]+'\t'+str(num[i])+'\t'+y+'\n')
	else:
		output.append(peak[i][0]+'\t'+peak[i][1]+'\t'+peak[i][2]+'\t'+peak[i][3]+'\t'+str(num[i])+'\t'+str(0)+'\n')

with open('reads.txt','w') as f:
	f.writelines(output)
f.close()
