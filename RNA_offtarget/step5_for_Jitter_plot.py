f1=open('merge_GATK_14sample_snp.txt','r')
f2=open('merge_GATK_10sample_Jitter.txt','w')
f2.write("Individual\tpercent\tStatus\tchr\tpos\tref_alt\n")

def refalt(gt1,gt2):
	change='NA'
	if gt1[0]==gt1[1] and gt2[0]==gt2[1]:
		# print([gt1,gt2])
		change=gt1[0]+'>'+gt2[0]
	if gt1[0]==gt1[1] and gt2[0] !=gt2[1]:
		#print([gt1,gt2])
		if gt1[0]==gt2[0]:
			change=gt1[0]+'>'+gt2[1]
		if gt1[0]==gt2[1]:
			change=gt1[0]+'>'+gt2[0]
	if gt1[0] !=gt1[1] and gt2[0] ==gt2[1]:
		# print([gt1,gt2])
		if gt2[0]==gt1[0]:
			change=gt1[1]+'>'+gt2[0]
		if gt2[0]==gt1[1]:
			change=gt1[0]+'>'+gt2[0]
		# print(change)
	if gt1[0] !=gt1[1] and gt2[0] !=gt2[1]:
		# print([gt1,gt2])
		if gt1[0]==gt2[0]:
			change=gt1[1]+'>'+gt2[1]
		if gt1[0]==gt2[1]:
			change=gt1[1]+'>'+gt2[0]
		if gt1[1]==gt2[0]:
			change=gt1[0]+'>'+gt2[1]
		if gt1[1]==gt2[1]:
			change=gt1[0]+'>'+gt2[0]
		# print(change)
	return change

dict1={}
same=[]
####headinfo="POS\tREF\t08431GFP\t08431neg\t08431pos\t08431WT\tNC1\tNC2\tNC3\tNC4\tNC5\tPC1\tPC2\tPC3\tPC4\tPC5\n"
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if len(set(tmp[2:6]))==1 and tmp[4] !='NA' and len(set(tmp[6:]))>1 and tmp[6:].count('NA')<1:
		if tmp[6] !=tmp[4] and tmp[6] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[6])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC1_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[7] !=tmp[4] and tmp[7] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[7])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC2_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[8] !=tmp[4] and tmp[8] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[8])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC3_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[9] !=tmp[4] and tmp[9] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[9])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC4_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[10] !=tmp[4] and tmp[10] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[10])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC5_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[11] !=tmp[4] and tmp[11] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[11])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC1_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[12] !=tmp[4] and tmp[12] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[12])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC2_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[13] !=tmp[4] and tmp[13] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[13])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC3_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[14] !=tmp[4] and tmp[14] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[14])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC4_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[15] !=tmp[4] and tmp[15] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[15])
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC5_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
	else:
		same.append(i)

for i in dict1.keys():
	list1=dict1[i]
	Status=i[5:7]
	Individual=i[5:8]
	f3=open('../'+i,'r')
	tmp_dict={}
	for j in f3.readlines()[1:]:
		tmp=j.strip().split('\t')
		pos='\t'.join([tmp[0],tmp[1]])
		chrm=tmp[0].split(':')[0]
		percent=str(round(float(tmp[4]), 4)*100)
		ref_alt=tmp[0]+'|'+tmp[1]+'>'+tmp[2]##2:96796421|A>G
		tmp_dict[pos]=[Individual,percent,Status,chrm,tmp[0],ref_alt]
	for j in list1:
		f2.write('\t'.join(tmp_dict[j])+'\n')
	f3.close()
f1.close()
f2.close()
#f4.close()
