# from collections import Counter
import numpy as np
f0=open('merge_26sample_filter_add_4method.txt','w')
f0a=open('merge_26sample_filter_add_4method_overlap.txt','w')
f0b=open('merge_26sample_filter_add_4method_overlap_venn.txt','w')
f0c=open('merge_26sample_filter_add_4method_overlapwithsamplename.txt','w')
f0d=open('merge_26sample_filter_add_4method_overlapwithsamplenamewithNC.txt','w')

lofreq={}
f1=open('merge_26sample_refgenome_adddonor_change_lofreq2.txt','r')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	str2='|'.join(tmp[0:2])
	lofreq.setdefault(tmp[2],[]).append(str1)
	f0b.write(tmp[2]+'\t'+str2+'\tlofreq2\n')
f1.close()

strelka={}
f1=open('merge_26sample_refgenome_adddonor_change_strelka2.txt','r')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	str2='|'.join(tmp[0:2])
	strelka.setdefault(tmp[2],[]).append(str1)
	f0b.write(tmp[2]+'\t'+str2+'\tstrelka\n')
f1.close()
platypus={}
f1=open('merge_26sample_refgenome_adddonor_change_platypus.txt','r')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	str2='|'.join(tmp[0:2])
	platypus.setdefault(tmp[2],[]).append(str1)
	f0b.write(tmp[2]+'\t'+str2+'\tplatypus\n')
f1.close()
gatk={}
f1=open('merge_26sample_refgenome_adddonor_change_gatk.txt','r')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	str2='|'.join(tmp[0:2])
	gatk.setdefault(tmp[2],[]).append(str1)
	f0b.write(tmp[2]+'\t'+str2+'\tgatk\n')
f1.close()

# 090202body	090202head	NC1	NC2	NC3	NC4	NC5	NC6	NC7	NC8	NC9	cloningJR	cloningPF	cloningW	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8
NC_sample=['090202body','090202head','NC1', 'NC2', 'NC3','NC4', 'NC5', 'NC6','NC7', 'NC8', 'NC9']
PC_sample=['cloningJR','cloningPF','cloningW','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8']

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
