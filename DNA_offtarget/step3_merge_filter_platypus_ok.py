f1=open('merge_15sample_add_platypus_v1.txt','r')
f2=open('merge_15sample_filter_platypus_v1.txt','r')
f3=open('merge_15sample_filter_add_platypus_v1.txt','w')
f4=open('merge_15sample_filter_add_count_platypus_v1.txt','w')
f5=open('merge_15sample_filter_add_change_platypus_v1.txt','w')

def compare(a,b):
	mycount=0
	for i in b:
		if i != a and i !='NA':
			mycount=mycount+1
	return mycount
filter_list=[]
filter_list1=[]
filter_list2=[]
same=[]
for i in f1.readlines():
	tmp=i.strip().split('\t')
	pc=compare(tmp[2],tmp[6:])
	nc=compare(tmp[2],tmp[3:6])
	if len(set(tmp[3:]))==1 or tmp[2] =='NA':#### or 'NA' in tmp:
		same.append(i)
	else:
		if pc>=0:####no filter
			filter_list1.append(tmp[0])
			filter_list.append(i.strip())
			f3.write(i)

f1.close()
for i in f2.readlines()[1:]:
	tmp=i.strip().split('\t')
	pc=compare(tmp[2],tmp[6:])
	nc=compare(tmp[2],tmp[3:6])
	if len(set(tmp[3:]))==1 or tmp[2] =='NA':#### or 'NA' in tmp:
		same.append(i)
	else:
		if pc>=0:####no filter
			filter_list2.append(tmp[0])
			filter_list.append(i.strip())
			f3.write(i)

print ('platypus',len(filter_list1))
print ('platypus',len(filter_list2))
f2.close()
difset1=list(set(filter_list2)-set(filter_list1))
difset2=list(set(filter_list1)-set(filter_list2))
print ('platypus',len(difset1))
print ('platypus',len(difset2))
print('platypus',len(set(filter_list)))
same=[]

f4.write("Individual\tbase_change\tnumber\n")
dict1={}

for i in list(set(filter_list)):
	tmp=i.strip().split('\t')
	if tmp[3] !=tmp[2] and tmp[3] !='NA' and tmp[2] !='NA' and (tmp[3] not in tmp[6:]):
		str1=tmp[2][1]+'>'+tmp[3][1]
		if tmp[2][1] ==tmp[3][1]:
			str1=tmp[2][0]+'>'+tmp[3][0]
			# print (i)
		dict1.setdefault('NC1',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[3],'NC1'])+'\n')
	if tmp[4] !=tmp[2] and tmp[4] !='NA' and tmp[2] !='NA' and (tmp[4] not in tmp[6:]):
		str1=tmp[2][1]+'>'+tmp[4][1]
		if tmp[2][1] ==tmp[4][1]:
			str1=tmp[2][0]+'>'+tmp[4][0]
		dict1.setdefault('NC2',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[4],'NC2'])+'\n')
	if tmp[5] !=tmp[2] and tmp[5] !='NA' and tmp[2] !='NA' and (tmp[5] not in tmp[6:]):
		str1=tmp[2][1]+'>'+tmp[5][1]
		if tmp[2][1] ==tmp[5][1]:
			str1=tmp[2][0]+'>'+tmp[5][0]
		dict1.setdefault('NC3',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[5],'NC3'])+'\n')
	if tmp[6] !=tmp[2] and tmp[6] !='NA' and tmp[2] !='NA' and (tmp[6] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[6][1]
		if tmp[2][1] ==tmp[6][1]:
			str1=tmp[2][0]+'>'+tmp[6][0]
		dict1.setdefault('PC1',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[6],'PC1'])+'\n')
	if tmp[7] !=tmp[2] and tmp[7] !='NA' and tmp[2] !='NA' and (tmp[7] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[7][1]
		if tmp[2][1] ==tmp[7][1]:
			str1=tmp[2][0]+'>'+tmp[7][0]
		dict1.setdefault('PC2',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[7],'PC2'])+'\n')
	if tmp[8] !=tmp[2] and tmp[8] !='NA' and tmp[2] !='NA' and (tmp[8] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[8][1]
		if tmp[2][1] ==tmp[8][1]:
			str1=tmp[2][0]+'>'+tmp[8][0]
		dict1.setdefault('PC3',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[8],'PC3'])+'\n')
	if tmp[9] !=tmp[2] and tmp[9] !='NA' and tmp[2] !='NA' and (tmp[9] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[9][1]
		if tmp[2][1] ==tmp[9][1]:
			str1=tmp[2][0]+'>'+tmp[9][0]
		dict1.setdefault('PC4',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[9],'PC4'])+'\n')
	if tmp[10] !=tmp[2] and tmp[10] !='NA' and tmp[2] !='NA' and (tmp[10] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[10][1]
		if tmp[2][1] ==tmp[10][1]:
			str1=tmp[2][0]+'>'+tmp[10][0]
		dict1.setdefault('PC5',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[10],'PC5'])+'\n')
	if tmp[11] !=tmp[2] and tmp[11] !='NA' and tmp[2] !='NA' and (tmp[11] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[11][1]
		if tmp[2][1] ==tmp[11][1]:
			str1=tmp[2][0]+'>'+tmp[11][0]
		dict1.setdefault('PC6',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[11],'PC6'])+'\n')
	if tmp[12] !=tmp[2] and tmp[12] !='NA' and tmp[2] !='NA' and (tmp[12] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[12][1]
		if tmp[2][1] ==tmp[12][1]:
			str1=tmp[2][0]+'>'+tmp[12][0]
		dict1.setdefault('PC7',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[12],'PC7'])+'\n')
	if tmp[13] !=tmp[2] and tmp[13] !='NA' and tmp[2] !='NA' and (tmp[13] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[13][1]
		if tmp[2][1] ==tmp[13][1]:
			str1=tmp[2][0]+'>'+tmp[13][0]
		dict1.setdefault('PC8',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[13],'PC8'])+'\n')
	if tmp[14] !=tmp[2] and tmp[14] !='NA' and tmp[2] !='NA' and (tmp[14] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[14][1]
		if tmp[2][1] ==tmp[14][1]:
			str1=tmp[2][0]+'>'+tmp[14][0]
		dict1.setdefault('cloningW',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[14],'cloningW'])+'\n')
	if tmp[15] !=tmp[2] and tmp[15] !='NA' and tmp[2] !='NA' and (tmp[15] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[15][1]
		if tmp[2][1] ==tmp[15][1]:
			str1=tmp[2][0]+'>'+tmp[15][0]
		dict1.setdefault('cloningPF',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[15],'cloningPF'])+'\n')
	if tmp[16] !=tmp[2] and tmp[16] !='NA' and tmp[2] !='NA' and (tmp[16] not in tmp[3:6]):
		str1=tmp[2][1]+'>'+tmp[16][1]
		if tmp[2][1] ==tmp[16][1]:
			str1=tmp[2][0]+'>'+tmp[16][0]
		dict1.setdefault('cloningJR',[]).append(str1)
		f5.write('\t'.join([tmp[0],str1,tmp[16],'cloningJR'])+'\n')

print (len(same))
for i in dict1.keys():
	list1=dict1[i]
	for j in list(set(list1)):
		f4.write(i+'\t'+j+'\t'+str(list1.count(j))+'\n')
f3.close()
f4.close()
f5.close()