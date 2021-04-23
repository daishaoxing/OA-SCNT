import os
import threading
from collections import Counter
import _pickle as cPickle
import numpy as np


dp=40
f1=open('07027_gatk_SNP.txt','r')
list07027=[]
dictref={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		list07027.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
f1.close()

f1=open('08431GFP_gatk_SNP.txt','r')
list08431GFP=[]
GFP_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		GFP_dict[tmp[0]]=tmp[3]
		list08431GFP.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		GFP_dict[tmp[0]]='NA'
f1.close()

f1=open('08431WT_gatk_SNP.txt','r')
list08431WT=[]
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		list08431WT.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
f1.close()

f1=open('NC1_gatk_SNP.txt','r')
listNC1=[]
NC1_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		NC1_dict[tmp[0]]=tmp[3]
		listNC1.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		NC1_dict[tmp[0]]='NA'
f1.close()

f1=open('NC2_gatk_SNP.txt','r')
listNC2=[]
NC2_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		NC2_dict[tmp[0]]=tmp[3]
		listNC2.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		NC2_dict[tmp[0]]='NA'
f1.close()

f1=open('NC3_gatk_SNP.txt','r')
listNC3=[]
NC3_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		NC3_dict[tmp[0]]=tmp[3]
		listNC3.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		NC3_dict[tmp[0]]='NA'
f1.close()

f1=open('PC1_gatk_SNP.txt','r')
listPC1=[]
PC1_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC1_dict[tmp[0]]=tmp[3]
		listPC1.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC1_dict[tmp[0]]='NA'
f1.close()

f1=open('PC2_gatk_SNP.txt','r')
listPC2=[]
PC2_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC2_dict[tmp[0]]=tmp[3]
		listPC2.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC2_dict[tmp[0]]='NA'
f1.close()

f1=open('PC3_gatk_SNP.txt','r')
listPC3=[]
PC3_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC3_dict[tmp[0]]=tmp[3]
		listPC3.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC3_dict[tmp[0]]='NA'
f1.close()

f1=open('PC4_gatk_SNP.txt','r')
listPC4=[]
PC4_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC4_dict[tmp[0]]=tmp[3]
		listPC4.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC4_dict[tmp[0]]='NA'
f1.close()

f1=open('PC5_gatk_SNP.txt','r')
listPC5=[]
PC5_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC5_dict[tmp[0]]=tmp[3]
		listPC5.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC5_dict[tmp[0]]='NA'
f1.close()

f1=open('PC6_gatk_SNP.txt','r')
listPC6=[]
PC6_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC6_dict[tmp[0]]=tmp[3]
		listPC6.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC6_dict[tmp[0]]='NA'
f1.close()

f1=open('PC7_gatk_SNP.txt','r')
listPC7=[]
PC7_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC7_dict[tmp[0]]=tmp[3]
		listPC7.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC7_dict[tmp[0]]='NA'
f1.close()

f1=open('PC8_gatk_SNP.txt','r')
listPC8=[]
PC8_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		PC8_dict[tmp[0]]=tmp[3]
		listPC8.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		PC8_dict[tmp[0]]='NA'
f1.close()

f1=open('cloning-W_gatk_SNP.txt','r')
listcloningW=[]
cloningW_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		cloningW_dict[tmp[0]]=tmp[3]
		listcloningW.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		cloningW_dict[tmp[0]]='NA'
f1.close()

f1=open('cloning-PF_gatk_SNP.txt','r')
listcloningPF=[]
cloningPF_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		cloningPF_dict[tmp[0]]=tmp[3]
		listcloningPF.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		cloningPF_dict[tmp[0]]='NA'
f1.close()

f1=open('cloning-JR_gatk_SNP.txt','r')
listcloningJR=[]
cloningJR_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		dictref[tmp[0]]=tmp[1]
		cloningJR_dict[tmp[0]]=tmp[3]
		listcloningJR.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
		cloningJR_dict[tmp[0]]='NA'
f1.close()
#############background 81 macaca DP>20
f1=open('/home/devdata/nyy/nyy_GFP/GFP/filter_vcf/population_fillterSNP_v2.txt','r')
population=[]
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	population.append(tmp[0])
f1.close()
background=list(set(population).union(list07027))

####same sites for all sample
overlap1=list(set(list08431GFP).intersection(listNC1,listNC2,listNC3,listcloningW,listcloningPF,listcloningJR,listPC1,listPC2,listPC3,listPC4,listPC5,listPC6,listPC7,listPC8))
####sites different with ref genome
overlap2=list(set(list08431GFP).intersection(list08431WT))
print ('gatk',population[8])##0
print ('gatk',listcloningJR[8])##0
print ('gatk',len(overlap1))##0
print ('gatk',len(overlap2))###3239058
difset1=set(overlap2)-set(overlap1)
difset2=list(set(difset1)-set(population))

print ('gatk',len(difset1))##3239058
print ('gatk',len(difset2))##3239058

fout=open('merge_15sample_add_gatk_v1.txt','w')

for i in difset2:
	pos=i.split('|')[0]
	ref=dictref[pos]
	GFP=ref+ref;NC1=ref+ref;NC2=ref+ref;NC3=ref+ref;PC1=ref+ref;PC2=ref+ref;PC3=ref+ref;
	PC4=ref+ref;PC5=ref+ref;PC6=ref+ref;PC7=ref+ref;PC8=ref+ref;
	cloningW=ref+ref;cloningPF=ref+ref;cloningJR=ref+ref
	if pos in GFP_dict.keys():
		GFP=GFP_dict[pos]
	if pos in NC1_dict.keys():
		NC1=NC1_dict[pos]
	if pos in NC2_dict.keys():
		NC2=NC2_dict[pos]
	if pos in NC3_dict.keys():
		NC3=NC3_dict[pos]
	if pos in PC1_dict.keys():
		PC1=PC1_dict[pos]
	if pos in PC2_dict.keys():
		PC2=PC2_dict[pos]
	if pos in PC3_dict.keys():
		PC3=PC3_dict[pos]
	if pos in PC4_dict.keys():
		PC4=PC4_dict[pos]
	if pos in PC5_dict.keys():
		PC5=PC5_dict[pos]
	if pos in PC6_dict.keys():
		PC6=PC6_dict[pos]
	if pos in PC7_dict.keys():
		PC7=PC7_dict[pos]
	if pos in PC8_dict.keys():
		PC8=PC8_dict[pos]
	if pos in cloningW_dict.keys():
		cloningW=cloningW_dict[pos]
	if pos in cloningPF_dict.keys():
		cloningPF=cloningPF_dict[pos]
	if pos in cloningJR_dict.keys():
		cloningJR=cloningJR_dict[pos]
	fout.write(pos+'\t'+ref+'\t'+GFP+'\t'+'\t'.join([NC1,NC2,NC3,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,cloningW,cloningPF,cloningJR])+'\n')
fout.close()
