import os
import _pickle as cPickle
fout=open('merge_26sample_filter_lofreq2_v1.txt','w')

fout.write("pos"+'\t'+"ref"+'\t'+'\t'.join(['08431GFP','08431neg','08431pos','08431WT','090202body','090202head','NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9','cloningJR','cloningPF','cloningW','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'])+'\n')

# POS     REF     ALT     GT      AF      AF2     DP      info
# 1:220   G       A       GA      0.6153846153846154      0.38461538461538464     78.0    30,48

dp=40
dictref={}
def write_get_info(fname):
	fh=open(fname,'r')
	ref={}
	dict_gtype={}
	alt=[]
	for i in fh.readlines()[1:]:
		tmp=i.strip().split('\t')
		if float(tmp[6])>=dp:
			str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
			ref[tmp[0]]=tmp[1]
			dict_gtype[tmp[0]]=tmp[3]
			alt.append(str1)
		else:
			ref[tmp[0]]=tmp[1]
			dict_gtype[tmp[0]]='NA'
	fh.close()
	return ref,dict_gtype,alt

ref_07027,dict_07027,list_07027=write_get_info('../vcf_2021/07027_lofreq2_SNP.txt')
ref_08431GFP,dict_08431GFP,list_08431GFP=write_get_info('../vcf_2021/08431GFP_lofreq2_SNP.txt')
ref_08431neg,dict_08431neg,list_08431neg=write_get_info('../vcf_2021/08431neg_lofreq2_SNP.txt')
ref_08431pos,dict_08431pos,list_08431pos=write_get_info('../vcf_2021/08431pos_lofreq2_SNP.txt')
ref_08431WT,dict_08431WT,list_08431WT=write_get_info('../vcf_2021/08431WT_lofreq2_SNP.txt')
ref_090202body,dict_090202body,list_090202body=write_get_info('../vcf_2021/090202body_lofreq2_SNP.txt')
ref_090202head,dict_090202head,list_090202head=write_get_info('../vcf_2021/090202head_lofreq2_SNP.txt')
ref_NC1,dict_NC1,list_NC1=write_get_info('../vcf_2021/NC1_lofreq2_SNP.txt')
ref_NC2,dict_NC2,list_NC2=write_get_info('../vcf_2021/NC2_lofreq2_SNP.txt')
ref_NC3,dict_NC3,list_NC3=write_get_info('../vcf_2021/NC3_lofreq2_SNP.txt')
ref_NC4,dict_NC4,list_NC4=write_get_info('../vcf_2021/NC4_lofreq2_SNP.txt')
ref_NC5,dict_NC5,list_NC5=write_get_info('../vcf_2021/NC5_lofreq2_SNP.txt')
ref_NC6,dict_NC6,list_NC6=write_get_info('../vcf_2021/NC6_lofreq2_SNP.txt')
ref_NC7,dict_NC7,list_NC7=write_get_info('../vcf_2021/NC7_lofreq2_SNP.txt')
ref_NC8,dict_NC8,list_NC8=write_get_info('../vcf_2021/NC8_lofreq2_SNP.txt')
ref_NC9,dict_NC9,list_NC9=write_get_info('../vcf_2021/NC9_lofreq2_SNP.txt')
ref_cloningJR,dict_cloningJR,list_cloningJR=write_get_info('../vcf_2021/cloning-JR_lofreq2_SNP.txt')
ref_cloningPF,dict_cloningPF,list_cloningPF=write_get_info('../vcf_2021/cloning-PF_lofreq2_SNP.txt')
ref_cloningW,dict_cloningW,list_cloningW=write_get_info('../vcf_2021/cloning-W_lofreq2_SNP.txt')
ref_PC1,dict_PC1,list_PC1=write_get_info('../vcf_2021/PC1_lofreq2_SNP.txt')
ref_PC2,dict_PC2,list_PC2=write_get_info('../vcf_2021/PC2_lofreq2_SNP.txt')
ref_PC3,dict_PC3,list_PC3=write_get_info('../vcf_2021/PC3_lofreq2_SNP.txt')
ref_PC4,dict_PC4,list_PC4=write_get_info('../vcf_2021/PC4_lofreq2_SNP.txt')
ref_PC5,dict_PC5,list_PC5=write_get_info('../vcf_2021/PC5_lofreq2_SNP.txt')
ref_PC6,dict_PC6,list_PC6=write_get_info('../vcf_2021/PC6_lofreq2_SNP.txt')
ref_PC7,dict_PC7,list_PC7=write_get_info('../vcf_2021/PC7_lofreq2_SNP.txt')
ref_PC8,dict_PC8,list_PC8=write_get_info('../vcf_2021/PC8_lofreq2_SNP.txt')

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

allSNP=list(set(list_090202body+list_090202head+list_NC1+list_NC2+list_NC3+list_NC4+list_NC5+list_NC6+list_NC7+list_NC8+list_NC9+list_PC1+list_PC2+list_PC3+list_PC4+list_PC5+list_PC6+list_PC7+list_PC8+list_cloningW+list_cloningPF+list_cloningJR))

allpos=[]
for i in allSNP:
	allpos.append(i.split('|')[0])
#cPickle.dump(depth40list,open("depth40_pos.pkl","wb"))
depth40list = cPickle.load(open("/nas01/nyy/vcffrom27/depth40_vcfpos.pkl","rb"))

print (len(allpos))
allpos1=list(set(allpos).intersection(depth40list))
print (len(allpos1))

for i in allpos1:
	pos=i
	ref=dictref[pos]
	GFP1=ref+ref;GFP2=ref+ref;GFP3=ref+ref;GFP4=ref+ref;
	NC1=ref+ref;NC2=ref+ref;NC3=ref+ref;NC4=ref+ref;NC5=ref+ref;NC6=ref+ref;NC7=ref+ref;NC8=ref+ref;NC9=ref+ref;S090202body=ref+ref;S090202head=ref+ref;
	PC1=ref+ref;PC2=ref+ref;PC3=ref+ref;PC4=ref+ref;PC5=ref+ref;PC6=ref+ref;PC7=ref+ref;PC8=ref+ref;cloningW=ref+ref;cloningPF=ref+ref;cloningJR=ref+ref
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
	fout.write(pos+'\t'+ref+'\t'+'\t'.join([GFP1,GFP2,GFP3,GFP4,S090202body,S090202head,NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,NC9,cloningJR,cloningPF,cloningW,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8])+'\n')
fout.close()

