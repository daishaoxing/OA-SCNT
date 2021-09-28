import os
import _pickle as cPickle
fout=open('merge_26sample_change_raw_gatk.txt','w')

fout.write("pso_base_change\tGT\tsample\n")
# POS     REF     ALT     GT      AF      AF2     DP      info
# 1:220   G       A       GA      0.6153846153846154      0.38461538461538464     78.0    30,48

dp=40
dictref={}
def write_get_info(fname):
	fh=open(fname,'r')
	dict1={}
	dict2={}
	for i in fh.readlines()[1:]:
		tmp=i.strip().split('\t')
		if float(tmp[6])>=dp:
			str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
			dict1[tmp[0]]=tmp[3]
			dict2[tmp[0]]=str1
	fh.close()
	return dict1,dict2

dict1_08431GFP,dict2_08431GFP=write_get_info('../vcf_2021/08431GFP_gatk_SNP.txt')
dict1_08431neg,dict2_08431neg=write_get_info('../vcf_2021/08431neg_gatk_SNP.txt')
dict1_08431pos,dict2_08431pos=write_get_info('../vcf_2021/08431pos_gatk_SNP.txt')
dict1_08431WT,dict2_08431WT=write_get_info('../vcf_2021/08431WT_gatk_SNP.txt')
dict1_090202body,dict2_090202body=write_get_info('../vcf_2021/090202body_gatk_SNP.txt')
dict1_090202head,dict2_090202head=write_get_info('../vcf_2021/090202head_gatk_SNP.txt')
dict1_NC1,dict2_NC1=write_get_info('../vcf_2021/NC1_gatk_SNP.txt')
dict1_NC2,dict2_NC2=write_get_info('../vcf_2021/NC2_gatk_SNP.txt')
dict1_NC3,dict2_NC3=write_get_info('../vcf_2021/NC3_gatk_SNP.txt')
dict1_NC4,dict2_NC4=write_get_info('../vcf_2021/NC4_gatk_SNP.txt')
dict1_NC5,dict2_NC5=write_get_info('../vcf_2021/NC5_gatk_SNP.txt')
dict1_NC6,dict2_NC6=write_get_info('../vcf_2021/NC6_gatk_SNP.txt')
dict1_NC7,dict2_NC7=write_get_info('../vcf_2021/NC7_gatk_SNP.txt')
dict1_NC8,dict2_NC8=write_get_info('../vcf_2021/NC8_gatk_SNP.txt')
dict1_NC9,dict2_NC9=write_get_info('../vcf_2021/NC9_gatk_SNP.txt')
dict1_cloningJR,dict2_cloningJR=write_get_info('../vcf_2021/cloning-JR_gatk_SNP.txt')
dict1_cloningPF,dict2_cloningPF=write_get_info('../vcf_2021/cloning-PF_gatk_SNP.txt')
dict1_cloningW,dict2_cloningW=write_get_info('../vcf_2021/cloning-W_gatk_SNP.txt')
dict1_PC1,dict2_PC1=write_get_info('../vcf_2021/PC1_gatk_SNP.txt')
dict1_PC2,dict2_PC2=write_get_info('../vcf_2021/PC2_gatk_SNP.txt')
dict1_PC3,dict2_PC3=write_get_info('../vcf_2021/PC3_gatk_SNP.txt')
dict1_PC4,dict2_PC4=write_get_info('../vcf_2021/PC4_gatk_SNP.txt')
dict1_PC5,dict2_PC5=write_get_info('../vcf_2021/PC5_gatk_SNP.txt')
dict1_PC6,dict2_PC6=write_get_info('../vcf_2021/PC6_gatk_SNP.txt')
dict1_PC7,dict2_PC7=write_get_info('../vcf_2021/PC7_gatk_SNP.txt')
dict1_PC8,dict2_PC8=write_get_info('../vcf_2021/PC8_gatk_SNP.txt')

#cPickle.dump(depth40list,open("depth40_pos.pkl","wb"))
depth40list = cPickle.load(open("/nas01/nyy/vcffrom27/depth40_vcfpos.pkl","rb"))

filter_pos=list(set(dict2_08431GFP.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_08431GFP[i],dict1_08431GFP[i]])+'\t08431GFP\n')

filter_pos=list(set(dict2_08431neg.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_08431neg[i],dict1_08431neg[i]])+'\t08431neg\n')

filter_pos=list(set(dict2_08431pos.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_08431pos[i],dict1_08431pos[i]])+'\t08431pos\n')

filter_pos=list(set(dict2_08431WT.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_08431WT[i],dict1_08431WT[i]])+'\t08431WT\n')

filter_pos=list(set(dict2_090202body.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_090202body[i],dict1_090202body[i]])+'\t090202body\n')

filter_pos=list(set(dict2_090202head.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_090202head[i],dict1_090202head[i]])+'\t090202head\n')

filter_pos=list(set(dict2_cloningJR.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_cloningJR[i],dict1_cloningJR[i]])+'\tcloningJR\n')

filter_pos=list(set(dict2_cloningPF.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_cloningPF[i],dict1_cloningPF[i]])+'\tcloningPF\n')

filter_pos=list(set(dict2_cloningW.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_cloningW[i],dict1_cloningW[i]])+'\tcloningW\n')

filter_pos=list(set(dict2_NC1.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC1[i],dict1_NC1[i]])+'\tNC1\n')

filter_pos=list(set(dict2_NC2.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC2[i],dict1_NC2[i]])+'\tNC2\n')

filter_pos=list(set(dict2_NC3.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC3[i],dict1_NC3[i]])+'\tNC3\n')

filter_pos=list(set(dict2_NC4.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC4[i],dict1_NC4[i]])+'\tNC4\n')

filter_pos=list(set(dict2_NC5.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC5[i],dict1_NC5[i]])+'\tNC5\n')

filter_pos=list(set(dict2_NC6.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC6[i],dict1_NC6[i]])+'\tNC6\n')

filter_pos=list(set(dict2_NC7.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC7[i],dict1_NC7[i]])+'\tNC7\n')

filter_pos=list(set(dict2_NC8.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC8[i],dict1_NC8[i]])+'\tNC8\n')

filter_pos=list(set(dict2_NC9.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_NC9[i],dict1_NC9[i]])+'\tNC9\n')

filter_pos=list(set(dict2_PC1.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC1[i],dict1_PC1[i]])+'\tPC1\n')

filter_pos=list(set(dict2_PC2.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC2[i],dict1_PC2[i]])+'\tPC2\n')

filter_pos=list(set(dict2_PC3.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC3[i],dict1_PC3[i]])+'\tPC3\n')

filter_pos=list(set(dict2_PC4.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC4[i],dict1_PC4[i]])+'\tPC4\n')

filter_pos=list(set(dict2_PC5.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC5[i],dict1_PC5[i]])+'\tPC5\n')

filter_pos=list(set(dict2_PC6.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC6[i],dict1_PC6[i]])+'\tPC6\n')

filter_pos=list(set(dict2_PC7.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC7[i],dict1_PC7[i]])+'\tPC7\n')

filter_pos=list(set(dict2_PC8.keys()).intersection(depth40list))
for i in filter_pos:
	fout.write('\t'.join([dict2_PC8[i],dict1_PC8[i]])+'\tPC8\n')

fout.close()
