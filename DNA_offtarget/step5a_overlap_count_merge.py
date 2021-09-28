f1=open('merge_26sample_filter_add_4method.txt','r')
f2=open('merge_26sample_filter_add_4method_count_merge.txt','w')
f2.write("Individual\tbase_change\tnumber\n")

def sumlist(a):
	b=0
	for i in a:
		b=b+int(i)
	return str(b)

dict1={}
for i in f1.readlines():
	tmp=i.strip().split('\t')
	if tmp[1]=='A>T' or tmp[1]=='T>A':
		str1=tmp[0]+'\t'+'A>T'
		dict1.setdefault(str1,[]).append(tmp[2])
	if tmp[1]=='T>C' or tmp[1]=='A>G':
		str1=tmp[0]+'\t'+'A>G'
		dict1.setdefault(str1,[]).append(tmp[2])
	if tmp[1]=='T>G' or tmp[1]=='A>C':
		str1=tmp[0]+'\t'+'A>C'
		dict1.setdefault(str1,[]).append(tmp[2])
	if tmp[1]=='C>A' or tmp[1]=='G>T':
		str1=tmp[0]+'\t'+'C>A'
		dict1.setdefault(str1,[]).append(tmp[2])
	if tmp[1]=='C>T' or tmp[1]=='G>A':
		str1=tmp[0]+'\t'+'C>T'
		dict1.setdefault(str1,[]).append(tmp[2])
	if tmp[1]=='C>G' or tmp[1]=='G>C':
		str1=tmp[0]+'\t'+'C>G'
		dict1.setdefault(str1,[]).append(tmp[2])
for i in dict1.keys():
	f2.write(i+'\t'+sumlist(dict1[i])+'\n')
f1.close()
f2.close()
