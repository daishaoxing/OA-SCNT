f1=open('merge_GATK_14sample_snp.txt','r')
f2=open('GATK_10sample_edit_site.txt','w')
dict1={}
same=[]


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

####headinfo="POS\tREF\t08431GFP\t08431neg\t08431pos\t08431WT\tNC1\tNC2\tNC3\tNC4\tNC5\tPC1\tPC2\tPC3\tPC4\tPC5\n"
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	# m=re.search(r'\|',i)
	refgene=tmp[1]+tmp[1]
	#if refgene==tmp[2]:
	if len(set(tmp[2:6]))==1 and tmp[4] !='NA' and len(set(tmp[6:]))>1 and tmp[6:].count('NA')<1:
		if tmp[6] !=tmp[4] and tmp[6] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[6])
			if str1=='T>C':
				dict1.setdefault('NC1',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC1',[]).append('|'.join([tmp[0],'+']))
		if tmp[7] !=tmp[4] and tmp[7] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[7])
			if str1=='T>C':
				dict1.setdefault('NC2',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC2',[]).append('|'.join([tmp[0],'+']))
		if tmp[8] !=tmp[4] and tmp[8] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[8])
			if str1=='T>C':
				dict1.setdefault('NC3',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC3',[]).append('|'.join([tmp[0],'+']))
		if tmp[9] !=tmp[4] and tmp[9] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[9])
			if str1=='T>C':
				dict1.setdefault('NC4',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC4',[]).append('|'.join([tmp[0],'+']))
		if tmp[10] !=tmp[4] and tmp[10] not in tmp[11:]:
			str1=refalt(tmp[4],tmp[10])
			if str1=='T>C':
				dict1.setdefault('NC5',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('NC5',[]).append('|'.join([tmp[0],'+']))
		if tmp[11] !=tmp[4] and tmp[11] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[11])
			if str1=='T>C':
				dict1.setdefault('PC1',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC1',[]).append('|'.join([tmp[0],'+']))
		if tmp[12] !=tmp[4] and tmp[12] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[12])
			if str1=='T>C':
				dict1.setdefault('PC2',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC2',[]).append('|'.join([tmp[0],'+']))
		if tmp[13] !=tmp[4] and tmp[13] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[13])
			if str1=='T>C':
				dict1.setdefault('PC3',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC3',[]).append('|'.join([tmp[0],'+']))
		if tmp[14] !=tmp[4] and tmp[14] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[14])
			if str1=='T>C':
				dict1.setdefault('PC4',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC4',[]).append('|'.join([tmp[0],'+']))
		if tmp[15] !=tmp[4] and tmp[15] not in tmp[6:11]:
			str1=refalt(tmp[4],tmp[15])
			if str1=='T>C':
				dict1.setdefault('PC5',[]).append('|'.join([tmp[0],'-']))
			if str1=='A>G':
				dict1.setdefault('PC5',[]).append('|'.join([tmp[0],'+']))
	else:
		same.append(i)
print (len(same))
for i in dict1.keys():
	list1=dict1[i]
	for j in list1:
		f2.write(i+'\t'+j+'\n')
f1.close()
f2.close()
