f1=open('merge_GATK_10sample_snp.txt','r')
f2=open('merge_GATK_10sample_Jitter_v1.txt','w')
f2.write("Individual\tpercent\tStatus\tchr\n")
dict1={}
same=[]
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	refgene=tmp[1]+tmp[1]
	if tmp[2] !='NA' and len(set(tmp[3:]))>1 and tmp[3:].count('NA')<1:##### and refgene==tmp[2]:
		if tmp[3] !='NA' and tmp[3] !=tmp[2] and tmp[3] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[3][1]
			if tmp[2][1] ==tmp[3][1]:
				str1=tmp[2][0]+'>'+tmp[3][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC1_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[4] !='NA' and tmp[4] !=tmp[2] and tmp[4] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[4][1]
			if tmp[2][1] ==tmp[4][1]:
				str1=tmp[2][0]+'>'+tmp[4][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC2_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[5] !='NA' and tmp[5] !=tmp[2] and tmp[5] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[5][1]
			if tmp[2][1] ==tmp[5][1]:
				str1=tmp[2][0]+'>'+tmp[5][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC3_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[6] !='NA' and tmp[6] !=tmp[2] and tmp[6] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[6][1]
			if tmp[2][1] ==tmp[6][1]:
				str1=tmp[2][0]+'>'+tmp[6][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC4_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[7] !='NA' and tmp[7] !=tmp[2] and tmp[7] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[7][1]
			if tmp[2][1] ==tmp[7][1]:
				str1=tmp[2][0]+'>'+tmp[7][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-NC5_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[8] !='NA' and tmp[8] !=tmp[2] and tmp[8] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[8][1]
			if tmp[2][1] ==tmp[8][1]:
				str1=tmp[2][0]+'>'+tmp[8][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC1_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[9] !='NA' and tmp[9] !=tmp[2] and tmp[9] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[9][1]
			if tmp[2][1] ==tmp[9][1]:
				str1=tmp[2][0]+'>'+tmp[9][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC2_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[10] !='NA' and tmp[10] !=tmp[2] and tmp[10] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[10][1]
			if tmp[2][1] ==tmp[10][1]:
				str1=tmp[2][0]+'>'+tmp[10][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC3_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[11] !='NA' and tmp[11] !=tmp[2] and tmp[11] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[11][1]
			if tmp[2][1] ==tmp[11][1]:
				str1=tmp[2][0]+'>'+tmp[11][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC4_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
		if tmp[12] !='NA' and tmp[12] !=tmp[2] and tmp[12] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[12][1]
			if tmp[2][1] ==tmp[12][1]:
				str1=tmp[2][0]+'>'+tmp[12][0]
			if str1=='T>C' or str1=='A>G':
				dict1.setdefault('SCNT-PC5_GATK_SNP05.txt',[]).append('\t'.join([tmp[0],tmp[1]]))
	else:
		same.append(i)
print (len(same))

# gfpref={}
# gfpref1={}
# f4=open('/home/devdata/nyy/nyy_GFP/GFP/08431GFP_strelka2_SNP05.txt','r')
# for i in f4.readlines()[1:]:
	# tmp=i.strip().split('\t')
	# gfpref[tmp[0]]=round(float(tmp[4]), 4)
	# gfpref1[tmp[0]]=tmp[3]

for i in dict1.keys():
	list1=dict1[i]
	Status=i[5:7]
	Individual=i[5:8]
	f3=open(i,'r')
	for j in f3.readlines()[1:]:
		tmp=j.strip().split('\t')
		pos='\t'.join([tmp[0],tmp[1]])
		chr=tmp[0].split(':')[0]
		percent=str(round(float(tmp[4]), 4)*100)
		if pos in list1:
			f2.write('\t'.join([Individual,percent,Status,chr])+'\n')
	f3.close()
f1.close()
f2.close()
#f4.close()
