import os
import threading
from collections import Counter
import _pickle as cPickle
import numpy as np

fout=open('merge_15sample_filter_platypus_v1.txt','w')
fout.write("pos"+'\t'+"change"+'\t'+"GFP"+'\t'+'\t'.join(['NC1','NC2','NC3','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','cloningW','cloningPF','cloningJR'])+'\n')

dp=40
f1=open('07027_platypus_SNP.txt','r')
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

f1=open('08431GFP_platypus_SNP.txt','r')
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

f1=open('08431WT_platypus_SNP.txt','r')
list08431WT=[]
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
		list08431WT.append(str1)
	else:
		dictref[tmp[0]]=tmp[1]
f1.close()

f1=open('NC1_platypus_SNP.txt','r')
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

f1=open('NC2_platypus_SNP.txt','r')
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

f1=open('NC3_platypus_SNP.txt','r')
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

f1=open('PC1_platypus_SNP.txt','r')
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

f1=open('PC2_platypus_SNP.txt','r')
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

f1=open('PC3_platypus_SNP.txt','r')
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

f1=open('PC4_platypus_SNP.txt','r')
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

f1=open('PC5_platypus_SNP.txt','r')
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

f1=open('PC6_platypus_SNP.txt','r')
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

f1=open('PC7_platypus_SNP.txt','r')
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

f1=open('PC8_platypus_SNP.txt','r')
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

f1=open('cloning-W_platypus_SNP.txt','r')
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

f1=open('cloning-PF_platypus_SNP.txt','r')
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

f1=open('cloning-JR_platypus_SNP.txt','r')
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

background=list(set(population).union(list07027,list08431WT,list08431GFP))
listNC1=list(set(listNC1)-set(background))
listNC2=list(set(listNC2)-set(background))
listNC3=list(set(listNC3)-set(background))
listPC1=list(set(listPC1)-set(background))
listPC2=list(set(listPC2)-set(background))
listPC3=list(set(listPC3)-set(background))
listPC4=list(set(listPC4)-set(background))
listPC5=list(set(listPC5)-set(background))
listPC6=list(set(listPC6)-set(background))
listPC7=list(set(listPC7)-set(background))
listPC8=list(set(listPC8)-set(background))
listcloningW=list(set(listcloningW)-set(background))
listcloningPF=list(set(listcloningPF)-set(background))
listcloningJR=list(set(listcloningJR)-set(background))

allSNP=list(set(listNC1+listNC2+listNC3+listPC1+listPC2+listPC3+listPC4+listPC5+listPC6+listPC7+listPC8+listcloningW+listcloningPF+listcloningJR))
print (len(allSNP))

for i in allSNP:
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



