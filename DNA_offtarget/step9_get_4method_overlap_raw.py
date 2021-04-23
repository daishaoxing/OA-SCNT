import os
fout=open('4method_overlap_raw_count_dp20_filter.txt','w')
list1=['07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR']
dp=20
#############background 81 macaca DP>20
f1=open('/home/devdata/nyy/nyy_GFP/GFP/filter_vcf/population_fillterSNP_v2.txt','r')
population=[]
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	population.append(tmp[0])
f1.close()

for i in list1:
	fname1=i+'_gatk_SNP.txt'
	fname2=i+'_platypus_SNP.txt'
	fname3=i+'_lofreq2_SNP.txt'
	fname4=i+'_strelka2_SNP.txt'
	gatk=[]
	platypus=[]
	lofreq=[]
	strelka=[]
	f1=open(fname1,'r',errors='ignore')
	for j in f1.readlines()[1:]:
		tmp=j.strip().split('\t')
		if float(tmp[6])>=dp:
			str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
			gatk.append(str1)
	f1.close()
	f1=open(fname2,'r',errors='ignore')
	for j in f1.readlines()[1:]:
		tmp=j.strip().split('\t')
		if float(tmp[6])>=dp:
			str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
			platypus.append(str1)
	f1.close()
	f1=open(fname3,'r',errors='ignore')
	for j in f1.readlines()[1:]:
		tmp=j.strip().split('\t')
		if float(tmp[6])>=dp:
			str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
			lofreq.append(str1)
	f1.close()
	f1=open(fname4,'r',errors='ignore')
	for j in f1.readlines()[1:]:
		tmp=j.strip().split('\t')
		if float(tmp[6])>=dp:
			str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
			strelka.append(str1)
	f1.close()
	overlap=list(set(lofreq).intersection(strelka,platypus,gatk))
	print(len(overlap))
	# print(len(population))
	overlap_filter=list(set(overlap)-set(population))
	fout.write(str(len(overlap_filter))+'\t'+i+'\n')
fout.close()