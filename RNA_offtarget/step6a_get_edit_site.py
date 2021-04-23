f1=open('merge_GATK_10sample_snp.txt','r')
f2=open('GATK_10sample_edit_sitev2.txt','w')
dict1={}
same=[]
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	if tmp[2] !='NA' and len(set(tmp[3:]))>1 and tmp[3:].count('NA')<1:
		if tmp[3] !='NA' and tmp[3] !=tmp[2] and tmp[3] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[3][1]
			if tmp[2][1] ==tmp[3][1]:
				str1=tmp[2][0]+'>'+tmp[3][0]
			if str1=='T>C':
				dict1.setdefault('NC1',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC1',[]).append('|'.join([tmp[0],'+']))
		if tmp[4] !='NA' and tmp[4] !=tmp[2] and tmp[4] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[4][1]
			if tmp[2][1] ==tmp[4][1]:
				str1=tmp[2][0]+'>'+tmp[4][0]
			if str1=='T>C':
				dict1.setdefault('NC2',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC2',[]).append('|'.join([tmp[0],'+']))
		if tmp[5] !='NA' and tmp[5] !=tmp[2] and tmp[5] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[5][1]
			if tmp[2][1] ==tmp[5][1]:
				str1=tmp[2][0]+'>'+tmp[5][0]
			if str1=='T>C':
				dict1.setdefault('NC3',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC3',[]).append('|'.join([tmp[0],'+']))
		if tmp[6] !='NA' and tmp[6] !=tmp[2] and tmp[6] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[6][1]
			if tmp[2][1] ==tmp[6][1]:
				str1=tmp[2][0]+'>'+tmp[6][0]
			if str1=='T>C':
				dict1.setdefault('NC4',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC4',[]).append('|'.join([tmp[0],'+']))
		if tmp[7] !='NA' and tmp[7] !=tmp[2] and tmp[7] not in tmp[8:]:
			str1=tmp[2][1]+'>'+tmp[7][1]
			if tmp[2][1] ==tmp[7][1]:
				str1=tmp[2][0]+'>'+tmp[7][0]
			if str1=='T>C':
				dict1.setdefault('NC5',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC5',[]).append('|'.join([tmp[0],'+']))
		if tmp[8] !='NA' and tmp[8] !=tmp[2] and tmp[8] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[8][1]
			if tmp[2][1] ==tmp[8][1]:
				str1=tmp[2][0]+'>'+tmp[8][0]
			if str1=='T>C':
				dict1.setdefault('PC1',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC1',[]).append('|'.join([tmp[0],'+']))
		if tmp[9] !='NA' and tmp[9] !=tmp[2] and tmp[9] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[9][1]
			if tmp[2][1] ==tmp[9][1]:
				str1=tmp[2][0]+'>'+tmp[9][0]
			if str1=='T>C':
				dict1.setdefault('PC2',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC2',[]).append('|'.join([tmp[0],'+']))
		if tmp[10] !='NA' and tmp[10] !=tmp[2] and tmp[10] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[10][1]
			if tmp[2][1] ==tmp[10][1]:
				str1=tmp[2][0]+'>'+tmp[10][0]
			if str1=='T>C':
				dict1.setdefault('PC3',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC3',[]).append('|'.join([tmp[0],'+']))
		if tmp[11] !='NA' and tmp[11] !=tmp[2] and tmp[11] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[11][1]
			if tmp[2][1] ==tmp[11][1]:
				str1=tmp[2][0]+'>'+tmp[11][0]
			if str1=='T>C':
				dict1.setdefault('PC4',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC4',[]).append('|'.join([tmp[0],'+']))
		if tmp[12] !='NA' and tmp[12] !=tmp[2] and tmp[12] not in tmp[3:8]:
			str1=tmp[2][1]+'>'+tmp[12][1]
			if tmp[2][1] ==tmp[12][1]:
				str1=tmp[2][0]+'>'+tmp[12][0]
			if str1=='T>C':
				dict1.setdefault('PC5',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC5',[]).append('|'.join([tmp[0],'+']))

# for i in dict1.keys():
	# list1=dict1[i]
	# f2.write(i+'\t'+str(len(list1))+'\t'+'\t'.join(list1)+'\n')

for i in dict1.keys():
	list1=dict1[i]
	for j in list1:
		f2.write(i+'\t'+j+'\n')
f1.close()
f2.close()
