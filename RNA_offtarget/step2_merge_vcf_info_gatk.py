import os
dict1={}
dict2={}
dict3={}
dict4={}
dict5={}
dict6={}
dict7={}
dict8={}
dict9={}
dict10={}
dict11={}
snplist=[]
f2=open('merge_GATK_10sample_snp.txt','w')
# list1=['08431GFP', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8']

headinfo="POS\tREF\t08431GFP\tNC1\tNC2\tNC3\tNC4\tNC5\tPC1\tPC2\tPC3\tPC4\tPC5\n"
f2.write(headinfo)

dp=40
f1=open('/home/devdata/nyy/nyy_GFP/GFP/08431GFP_strelka2_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t') ###"POS\tREF\tALT\tGT\tAF\tAF2\tDP\tinfo\n"
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict1[str1]=tmp[3]
	else:
		dict1[str1]='NA'
f1.close()

f1=open('SCNT-NC1_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t') ##"POS\tREF\tALT\tGT\tAF\tAF2\tDP\n"
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict2[str1]=tmp[3]
	else:
		dict2[str1]='NA'
f1.close()

f1=open('SCNT-NC2_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict3[str1]=tmp[3]
	else:
		dict3[str1]='NA'
f1.close()

f1=open('SCNT-NC3_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict4[str1]=tmp[3]
	else:
		dict4[str1]='NA'
f1.close()

f1=open('SCNT-NC4_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict5[str1]=tmp[3]
	else:
		dict5[str1]='NA'
f1.close()

f1=open('SCNT-NC5_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict6[str1]=tmp[3]
	else:
		dict6[str1]='NA'
f1.close()

f1=open('SCNT-PC1_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict7[str1]=tmp[3]
	else:
		dict7[str1]='NA'
f1.close()

f1=open('SCNT-PC2_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict8[str1]=tmp[3]
	else:
		dict8[str1]='NA'
f1.close()

f1=open('SCNT-PC3_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict9[str1]=tmp[3]
	else:
		dict9[str1]='NA'
f1.close()

f1=open('SCNT-PC4_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict10[str1]=tmp[3]
	else:
		dict10[str1]='NA'
f1.close()

f1=open('SCNT-PC5_GATK_SNP05.txt','r')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	str1='\t'.join(tmp[0:2])
	snplist.append(str1)
	if float(tmp[6])>=dp:## and float(tmp[4])>0.05
		dict11[str1]=tmp[3]
	else:
		dict11[str1]='NA'
f1.close()

for i in set(snplist):
	(a,b)=i.split('\t')
	GFP=b+b;NC1=b+b;NC2=b+b;NC3=b+b;NC4=b+b;NC5=b+b;
	PC1=b+b;PC2=b+b;PC3=b+b;PC4=b+b;PC5=b+b
	if i in dict1.keys():
		GFP=dict1[i]
	if i in dict2.keys():
		NC1=dict2[i]
	if i in dict3.keys():
		NC2=dict3[i]
	if i in dict4.keys():
		NC3=dict4[i]
	if i in dict5.keys():
		NC4=dict5[i]
	if i in dict6.keys():
		NC5=dict6[i]
	if i in dict7.keys():
		PC1=dict7[i]
	if i in dict8.keys():
		PC2=dict8[i]
	if i in dict9.keys():
		PC3=dict9[i]
	if i in dict10.keys():
		PC4=dict10[i]
	if i in dict11.keys():
		PC5=dict11[i]
	f2.write(i+'\t'+'\t'.join([GFP,NC1,NC2,NC3,NC4,NC5,PC1,PC2,PC3,PC4,PC5])+'\n')

f2.close()
