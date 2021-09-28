fin=open('merge_26sample_filter_strelka2_v1.txt','r')
fout1=open('merge_26sample_refgenome_change_strelka2.txt','w')
fout2=open('merge_26sample_adddonor_change_strelka2.txt','w')
fout1.write("pso\tbase_change\tGT\tsample\n")##4:30262373      A>T     AT      NC4
fout2.write("pso\tbase_change\tGT\tsample\n")
dict1={}

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

for i in fin.readlines()[1:]:
	tmp=i.strip().split('\t')
	refGT=tmp[1]+tmp[1]
	if len(set(tmp[2:6]))==1 and tmp[4] !='NA' and len(set(tmp[2:])) >1:
	# if tmp[4] !='NA' and len(set(tmp[2:])) >1:
		if tmp[4]==refGT:
			if tmp[6] !=tmp[4] and tmp[6] !='NA' and (tmp[6] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[6])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[6],'090202body'])+'\n')
			if tmp[7] !=tmp[4] and tmp[7] !='NA' and (tmp[7] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[7])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[7],'090202head'])+'\n')
			if tmp[8] !=tmp[4] and tmp[8] !='NA' and (tmp[8] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[8])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[8],'NC1'])+'\n')
			if tmp[9] !=tmp[4] and tmp[9] !='NA' and (tmp[9] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[9])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[9],'NC2'])+'\n')
			if tmp[10] !=tmp[4] and tmp[10] !='NA' and (tmp[10] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[10])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[10],'NC3'])+'\n')
			if tmp[11] !=tmp[4] and tmp[11] !='NA' and (tmp[11] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[11])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[11],'NC4'])+'\n')
			if tmp[12] !=tmp[4] and tmp[12] !='NA' and (tmp[12] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[12])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[12],'NC5'])+'\n')
			if tmp[13] !=tmp[4] and tmp[13] !='NA' and (tmp[13] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[13])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[13],'NC6'])+'\n')
			if tmp[14] !=tmp[4] and tmp[14] !='NA' and (tmp[14] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[14])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[14],'NC7'])+'\n')
			if tmp[15] !=tmp[4] and tmp[15] !='NA' and (tmp[15] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[15])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[15],'NC8'])+'\n')
			if tmp[16] !=tmp[4] and tmp[16] !='NA' and (tmp[16] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[16])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[16],'NC9'])+'\n')
##################SCNT-BE
			if tmp[17] !=tmp[4] and tmp[17] !='NA' and (tmp[17] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[17])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[17],'cloningJR'])+'\n')
			if tmp[18] !=tmp[4] and tmp[18] !='NA' and (tmp[18] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[18])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[18],'cloningPF'])+'\n')
			if tmp[19] !=tmp[4] and tmp[19] !='NA' and (tmp[19] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[19])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[19],'cloningW'])+'\n')
			if tmp[20] !=tmp[4] and tmp[20] !='NA' and (tmp[20] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[20])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[20],'PC1'])+'\n')
			if tmp[21] !=tmp[4] and tmp[21] !='NA' and (tmp[21] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[21])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[21],'PC2'])+'\n')
			if tmp[22] !=tmp[4] and tmp[22] !='NA' and (tmp[22] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[22])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[22],'PC3'])+'\n')
			if tmp[23] !=tmp[4] and tmp[23] !='NA' and (tmp[23] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[23])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[23],'PC4'])+'\n')
			if tmp[24] !=tmp[4] and tmp[24] !='NA' and (tmp[24] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[24])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[24],'PC5'])+'\n')
			if tmp[25] !=tmp[4] and tmp[25] !='NA' and (tmp[25] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[25])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[25],'PC6'])+'\n')
			if tmp[26] !=tmp[4] and tmp[26] !='NA' and (tmp[26] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[26])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[26],'PC7'])+'\n')
			if tmp[27] !=tmp[4] and tmp[27] !='NA' and (tmp[27] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[27])
				if str1 !='NA':
					fout1.write('\t'.join([tmp[0],str1,tmp[27],'PC8'])+'\n')
		else:
			if tmp[6] !=tmp[4] and tmp[6] !='NA' and (tmp[6] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[6])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[6],'090202body'])+'\n')
			if tmp[7] !=tmp[4] and tmp[7] !='NA' and (tmp[7] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[7])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[7],'090202head'])+'\n')
			if tmp[8] !=tmp[4] and tmp[8] !='NA' and (tmp[8] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[8])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[8],'NC1'])+'\n')
			if tmp[9] !=tmp[4] and tmp[9] !='NA' and (tmp[9] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[9])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[9],'NC2'])+'\n')
			if tmp[10] !=tmp[4] and tmp[10] !='NA' and (tmp[10] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[10])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[10],'NC3'])+'\n')
			if tmp[11] !=tmp[4] and tmp[11] !='NA' and (tmp[11] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[11])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[11],'NC4'])+'\n')
			if tmp[12] !=tmp[4] and tmp[12] !='NA' and (tmp[12] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[12])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[12],'NC5'])+'\n')
			if tmp[13] !=tmp[4] and tmp[13] !='NA' and (tmp[13] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[13])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[13],'NC6'])+'\n')
			if tmp[14] !=tmp[4] and tmp[14] !='NA' and (tmp[14] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[14])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[14],'NC7'])+'\n')
			if tmp[15] !=tmp[4] and tmp[15] !='NA' and (tmp[15] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[15])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[15],'NC8'])+'\n')
			if tmp[16] !=tmp[4] and tmp[16] !='NA' and (tmp[16] not in tmp[17:]):
				str1=refalt(tmp[4],tmp[16])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[16],'NC9'])+'\n')
##################SCNT-BE
			if tmp[17] !=tmp[4] and tmp[17] !='NA' and (tmp[17] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[17])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[17],'cloningJR'])+'\n')
			if tmp[18] !=tmp[4] and tmp[18] !='NA' and (tmp[18] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[18])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[18],'cloningPF'])+'\n')
			if tmp[19] !=tmp[4] and tmp[19] !='NA' and (tmp[19] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[19])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[19],'cloningW'])+'\n')
			if tmp[20] !=tmp[4] and tmp[20] !='NA' and (tmp[20] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[20])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[20],'PC1'])+'\n')
			if tmp[21] !=tmp[4] and tmp[21] !='NA' and (tmp[21] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[21])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[21],'PC2'])+'\n')
			if tmp[22] !=tmp[4] and tmp[22] !='NA' and (tmp[22] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[22])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[22],'PC3'])+'\n')
			if tmp[23] !=tmp[4] and tmp[23] !='NA' and (tmp[23] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[23])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[23],'PC4'])+'\n')
			if tmp[24] !=tmp[4] and tmp[24] !='NA' and (tmp[24] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[24])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[24],'PC5'])+'\n')
			if tmp[25] !=tmp[4] and tmp[25] !='NA' and (tmp[25] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[25])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[25],'PC6'])+'\n')
			if tmp[26] !=tmp[4] and tmp[26] !='NA' and (tmp[26] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[26])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[26],'PC7'])+'\n')
			if tmp[27] !=tmp[4] and tmp[27] !='NA' and (tmp[27] not in tmp[6:17]):
				str1=refalt(tmp[4],tmp[27])
				if str1 !='NA':
					fout2.write('\t'.join([tmp[0],str1,tmp[27],'PC8'])+'\n')
