fin=open('merge_26sample_filter_add_4method_overlapwithsamplename.txt','r')
# fin=open('merge_15sample_filter_add_4method_overlap.txt','r')
fout=open('4method_overlap_site_freq.txt','w')

dp=40
dictref={}

def write_get_info(fname):
	fh=open(fname,'r')
	ref={}
	dict_AF={}
	for i in fh.readlines()[1:]:
		tmp=i.strip().split('\t')
		if float(tmp[6])>=dp:
			str1=tmp[3]+'|'+str(round(float(tmp[4]),2))
			ref[tmp[0]]=tmp[1]
			dict_AF[tmp[0]]=str1
		else:
			ref[tmp[0]]=tmp[1]
	fh.close()
	return ref,dict_AF

ref_07027,dict_07027=write_get_info('../vcf_2021/07027_gatk_SNP.txt')
ref_08431GFP,dict_08431GFP=write_get_info('../vcf_2021/08431GFP_gatk_SNP.txt')
ref_08431neg,dict_08431neg=write_get_info('../vcf_2021/08431neg_gatk_SNP.txt')
ref_08431pos,dict_08431pos=write_get_info('../vcf_2021/08431pos_gatk_SNP.txt')
ref_08431WT,dict_08431WT=write_get_info('../vcf_2021/08431WT_gatk_SNP.txt')
ref_090202body,dict_090202body=write_get_info('../vcf_2021/090202body_gatk_SNP.txt')
ref_090202head,dict_090202head=write_get_info('../vcf_2021/090202head_gatk_SNP.txt')
ref_NC1,dict_NC1=write_get_info('../vcf_2021/NC1_gatk_SNP.txt')
ref_NC2,dict_NC2=write_get_info('../vcf_2021/NC2_gatk_SNP.txt')
ref_NC3,dict_NC3=write_get_info('../vcf_2021/NC3_gatk_SNP.txt')
ref_NC4,dict_NC4=write_get_info('../vcf_2021/NC4_gatk_SNP.txt')
ref_NC5,dict_NC5=write_get_info('../vcf_2021/NC5_gatk_SNP.txt')
ref_NC6,dict_NC6=write_get_info('../vcf_2021/NC6_gatk_SNP.txt')
ref_NC7,dict_NC7=write_get_info('../vcf_2021/NC7_gatk_SNP.txt')
ref_NC8,dict_NC8=write_get_info('../vcf_2021/NC8_gatk_SNP.txt')
ref_NC9,dict_NC9=write_get_info('../vcf_2021/NC9_gatk_SNP.txt')
ref_cloningJR,dict_cloningJR=write_get_info('../vcf_2021/cloning-JR_gatk_SNP.txt')
ref_cloningPF,dict_cloningPF=write_get_info('../vcf_2021/cloning-PF_gatk_SNP.txt')
ref_cloningW,dict_cloningW=write_get_info('../vcf_2021/cloning-W_gatk_SNP.txt')
ref_PC1,dict_PC1=write_get_info('../vcf_2021/PC1_gatk_SNP.txt')
ref_PC2,dict_PC2=write_get_info('../vcf_2021/PC2_gatk_SNP.txt')
ref_PC3,dict_PC3=write_get_info('../vcf_2021/PC3_gatk_SNP.txt')
ref_PC4,dict_PC4=write_get_info('../vcf_2021/PC4_gatk_SNP.txt')
ref_PC5,dict_PC5=write_get_info('../vcf_2021/PC5_gatk_SNP.txt')
ref_PC6,dict_PC6=write_get_info('../vcf_2021/PC6_gatk_SNP.txt')
ref_PC7,dict_PC7=write_get_info('../vcf_2021/PC7_gatk_SNP.txt')
ref_PC8,dict_PC8=write_get_info('../vcf_2021/PC8_gatk_SNP.txt')

dictref=ref_07027
dictref.update(ref_08431GFP)
dictref.update(ref_08431neg)
dictref.update(ref_08431pos)
dictref.update(ref_08431WT)
dictref.update(ref_090202body)
dictref.update(ref_090202head)
dictref.update(ref_NC1)
dictref.update(ref_NC2)
dictref.update(ref_NC3)
dictref.update(ref_NC4)
dictref.update(ref_NC5)
dictref.update(ref_NC6)
dictref.update(ref_NC7)
dictref.update(ref_NC8)
dictref.update(ref_NC9)
dictref.update(ref_cloningJR)
dictref.update(ref_cloningPF)
dictref.update(ref_cloningW)
dictref.update(ref_PC1)
dictref.update(ref_PC2)
dictref.update(ref_PC3)
dictref.update(ref_PC4)
dictref.update(ref_PC5)
dictref.update(ref_PC6)
dictref.update(ref_PC7)
dictref.update(ref_PC8)

fout.write('|'.join(['position','change','GT','num'])+'\t'+'refbase'+'\t'+'\t'.join(['08431GFP','08431neg','08431pos','08431WT','090202body','090202head','NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9','cloningJR','cloningPF','cloningW','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'])+'\n')

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
	GFP1=ref_frq;GFP2=ref_frq;GFP3=ref_frq;GFP4=ref_frq;
	NC1=ref_frq;NC2=ref_frq;NC3=ref_frq;NC4=ref_frq;NC5=ref_frq;NC6=ref_frq;NC7=ref_frq;NC8=ref_frq;NC9=ref_frq;S090202body=ref_frq;S090202head=ref_frq;
	PC1=ref_frq;PC2=ref_frq;PC3=ref_frq;PC4=ref_frq;PC5=ref_frq;PC6=ref_frq;PC7=ref_frq;PC8=ref_frq;cloningW=ref_frq;cloningPF=ref_frq;cloningJR=ref_frq
	if pos in dict_08431GFP.keys():
		GFP1=dict_08431GFP[pos]
	if pos in dict_08431neg.keys():
		GFP2=dict_08431neg[pos]
	if pos in dict_08431pos.keys():
		GFP3=dict_08431pos[pos]
	if pos in dict_08431WT.keys():
		GFP4=dict_08431WT[pos]
##########SCNT
	if pos in dict_090202body.keys():
		S090202body=dict_090202body[pos]
	if pos in dict_090202head.keys():
		S090202head=dict_090202head[pos]
	if pos in dict_NC1.keys():
		NC1=dict_NC1[pos]
	if pos in dict_NC2.keys():
		NC2=dict_NC2[pos]
	if pos in dict_NC3.keys():
		NC3=dict_NC3[pos]
	if pos in dict_NC4.keys():
		NC4=dict_NC4[pos]
	if pos in dict_NC5.keys():
		NC5=dict_NC5[pos]
	if pos in dict_NC6.keys():
		NC6=dict_NC6[pos]
	if pos in dict_NC7.keys():
		NC7=dict_NC7[pos]
	if pos in dict_NC8.keys():
		NC8=dict_NC8[pos]
	if pos in dict_NC9.keys():
		NC9=dict_NC9[pos]
##########SCNT-BE
	if pos in dict_cloningJR.keys():
		cloningJR=dict_cloningJR[pos]
	if pos in dict_cloningPF.keys():
		cloningPF=dict_cloningPF[pos]
	if pos in dict_cloningW.keys():
		cloningW=dict_cloningW[pos]
	if pos in dict_PC1.keys():
		PC1=dict_PC1[pos]
	if pos in dict_PC2.keys():
		PC2=dict_PC2[pos]
	if pos in dict_PC3.keys():
		PC3=dict_PC3[pos]
	if pos in dict_PC4.keys():
		PC4=dict_PC4[pos]
	if pos in dict_PC5.keys():
		PC5=dict_PC5[pos]
	if pos in dict_PC6.keys():
		PC6=dict_PC6[pos]
	if pos in dict_PC7.keys():
		PC7=dict_PC7[pos]
	if pos in dict_PC8.keys():
		PC8=dict_PC8[pos]
	fout.write('|'.join(tmp)+'|'+'|'.join(dict_sample[i])+'\t'+ref+'\t'+'\t'.join([GFP1,GFP2,GFP3,GFP4,S090202body,S090202head,NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,NC9,cloningJR,cloningPF,cloningW,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8])+'\n')

fin.close()
fout.close()


