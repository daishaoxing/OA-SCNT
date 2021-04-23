# from collections import Counter
import numpy as np
f0=open('merge_15sample_filter_add_4method.txt','w')
f0a=open('merge_15sample_filter_add_4method_overlap.txt','w')
f0b=open('merge_15sample_filter_add_4method_overlap_venn.txt','w')
f0c=open('merge_15sample_filter_add_4method_overlapwithsamplename.txt','w')
f0d=open('merge_15sample_filter_add_4method_overlapwithsamplenamewithNC.txt','w')

lofreq={}
f1=open('merge_15sample_filter_add_change_lofreq2_v1.txt','r')
for j in f1.readlines():
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:3])
	str2='|'.join(tmp[0:3])
	lofreq.setdefault(tmp[3],[]).append(str1)
	f0b.write(tmp[3]+'\t'+str2+'\tlofreq2\n')
f1.close()

strelka={}
f1=open('merge_15sample_filter_add_change_strelka2_v1.txt','r')
for j in f1.readlines():
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:3])
	str2='|'.join(tmp[0:3])
	strelka.setdefault(tmp[3],[]).append(str1)
	f0b.write(tmp[3]+'\t'+str2+'\tstrelka\n')
f1.close()
platypus={}
f1=open('merge_15sample_filter_add_change_platypus_v1.txt','r')
for j in f1.readlines():
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:3])
	str2='|'.join(tmp[0:3])
	platypus.setdefault(tmp[3],[]).append(str1)
	f0b.write(tmp[3]+'\t'+str2+'\tplatypus\n')
f1.close()
gatk={}
f1=open('merge_15sample_filter_add_change_gatk_v1.txt','r')
for j in f1.readlines():
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:3])
	str2='|'.join(tmp[0:3])
	gatk.setdefault(tmp[3],[]).append(str1)
	f0b.write(tmp[3]+'\t'+str2+'\tgatk\n')
f1.close()
NC_sample=['NC1', 'NC2', 'NC3']
PC_sample=['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloningW','cloningPF','cloningJR']

dict1={}
PCoverlap=[]
NCoverlap=[]
NCall=[]
for i in NC_sample:
	overlap=list(set(lofreq[i]).intersection(strelka[i],platypus[i],gatk[i]))
	all=list(set(lofreq[i]).union(strelka[i],platypus[i],gatk[i]))
	print (i,len(overlap))
	for j in all:
		NCall.append(j)
	for j in overlap:
		tmp=j.strip().split('\t')
		NCoverlap.append(j)
		f0d.write(j+'\t'+i+'\n')
		dict1.setdefault(i,[]).append(tmp[1])
for i in PC_sample:
	overlap=list(set(lofreq[i]).intersection(strelka[i],platypus[i],gatk[i]))
	# overlap=list(set(lofreq[i]).intersection(strelka[i],gatk[i]))
	print (i,len(overlap))
	for j in overlap:
		tmp=j.strip().split('\t')
		PCoverlap.append(j)
		f0c.write(j+'\t'+i+'\n')
		f0d.write(j+'\t'+i+'\n')
		dict1.setdefault(i,[]).append(tmp[1])

# PC18_count=dict(Counter(PCoverlap))
# PC18_count=sorted(PC18_count.items(), key=lambda d: d[1],reverse=True)

def unique_count(lst):
	lst = np.array(lst)
	(a,b)=np.unique(lst, return_counts=True)
	# return dict(zip(a,b))
	return list(zip(a,b))
PC18_count = unique_count(PCoverlap)
print(PC18_count[1])
# for i in PC18_count:
	# if i[0] not in NCall:
		# f0a.write(i[0]+'\t'+str(i[1])+'\n')
for i in PC18_count:
	f0a.write(i[0]+'\t'+str(i[1])+'\n')

for i in dict1.keys():
	list1=dict1[i]
	print(i,len(list1))
	for j in list(set(list1)):
		f0.write(i+'\t'+j+'\t'+str(list1.count(j))+'\n')
f0.close()
f0a.close()
f0b.close()
