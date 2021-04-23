fin=open('merge_15sample_filter_add_4method_overlapwithsamplename.txt','r')
# fin=open('merge_15sample_filter_add_4method_overlap.txt','r')
fout=open('4method_overlap_site_freq.txt','w')

dp=20
dictref={}
f1=open('08431GFP_gatk_SNP.txt','r')
list08431GFP=[]
GFP_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		GFP_dict[tmp[0]]=str1
f1.close()

f1=open('NC1_gatk_SNP.txt','r')
NC1_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		NC1_dict[tmp[0]]=str1
f1.close()

f1=open('NC2_gatk_SNP.txt','r')
NC2_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		NC2_dict[tmp[0]]=str1
f1.close()

f1=open('NC3_gatk_SNP.txt','r')
NC3_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		NC3_dict[tmp[0]]=str1
f1.close()

f1=open('PC1_gatk_SNP.txt','r')

PC1_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC1_dict[tmp[0]]=str1
f1.close()

f1=open('PC2_gatk_SNP.txt','r')
PC2_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC2_dict[tmp[0]]=str1
f1.close()

f1=open('PC3_gatk_SNP.txt','r')
PC3_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC3_dict[tmp[0]]=str1
f1.close()

f1=open('PC4_gatk_SNP.txt','r')
PC4_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC4_dict[tmp[0]]=str1
f1.close()
f1=open('PC5_gatk_SNP.txt','r')
PC5_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC5_dict[tmp[0]]=str1
f1.close()

f1=open('PC6_gatk_SNP.txt','r')
PC6_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC6_dict[tmp[0]]=str1
f1.close()

f1=open('PC7_gatk_SNP.txt','r')
PC7_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC7_dict[tmp[0]]=str1
f1.close()

f1=open('PC8_gatk_SNP.txt','r')
PC8_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		PC8_dict[tmp[0]]=str1
f1.close()

f1=open('cloning-W_gatk_SNP.txt','r')
cloningW_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		cloningW_dict[tmp[0]]=str1
f1.close()

f1=open('cloning-PF_gatk_SNP.txt','r')
cloningPF_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		cloningPF_dict[tmp[0]]=str1
f1.close()

f1=open('cloning-JR_gatk_SNP.txt','r')
cloningJR_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if float(tmp[6])>=dp:
		str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
		dictref[tmp[0]]=tmp[1]
		cloningJR_dict[tmp[0]]=str1
f1.close()

fout.write('|'.join(['position','change','GT','num'])+'\t'+'refbase'+'\t'+'GFP'+'\t'+'\t'.join(['NC1','NC2','NC3','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','cloningW','cloningPF','cloningJR'])+'\n')

lines=fin.readlines()
dict_sample={}
for i in lines:
	tmp=i.strip().split('\t')
	posinfo='\t'.join(tmp[0:-1])
	dict_sample.setdefault(posinfo,[]).append(tmp[-1])

for i in dict_sample.keys():
	tmp=i.strip().split('\t')
	pos=tmp[0]
	ref=dictref[pos]
	ref_frq=ref+ref+'|1'
	# GFP='NA';NC1='NA';NC2='NA';NC3='NA';PC1='NA';PC2='NA';PC3='NA';
	# PC4='NA';PC5='NA';PC6='NA';PC7='NA';PC8='NA';
	# cloningW='NA';cloningPF='NA';cloningJR='NA'
	GFP=ref_frq;NC1=ref_frq;NC2=ref_frq;NC3=ref_frq;PC1=ref_frq;PC2=ref_frq;PC3=ref_frq;
	PC4=ref_frq;PC5=ref_frq;PC6=ref_frq;PC7=ref_frq;PC8=ref_frq;
	cloningW=ref_frq;cloningPF=ref_frq;cloningJR=ref_frq
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
	fout.write('|'.join(tmp)+'|'+'|'.join(dict_sample[i])+'\t'+ref+'\t'+GFP+'\t'+'\t'.join([NC1,NC2,NC3,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,cloningW,cloningPF,cloningJR])+'\n')
fin.close()
fout.close()


