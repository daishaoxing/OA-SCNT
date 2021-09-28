f1=open('merge_26sample_filter_add_4method_count_merge.txt','r')
f2=open('merge_26sample_filter_add_4method_count_merge_proportion.txt','w')
f2.write("Individual\tbase_change\tnumber\tproportion\n")

lines=f1.readlines()[1:]
f1.close()
types=[]
samples=[]
sampletype=[]
dict1={}
dict1_sum={}
for i in lines:
	tmp=i.strip().split('\t')
	types.append(tmp[1])
	samples.append(tmp[0])
	dict1.setdefault(tmp[0],[]).append(int(tmp[2]))
types=list(set(types))
samples=list(set(samples))
for i in samples:
	dict1_sum[i]=sum(dict1[i])
	for j in types:
		sampletype.append(i+'_'+j)
sampletype2=[]
for i in lines:
	tmp=i.strip().split('\t')
	pro=float(tmp[-1])/dict1_sum[tmp[0]]
	sampletype2.append(tmp[0]+'_'+tmp[1])
	f2.write(tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+str(round(pro,4))+'\n')
for i in sampletype:
	if i not in sampletype2:
		(a,b)=i.split('_')
		f2.write(a+'\t'+b+'\t0\t0'+'\n')
f2.close()

####
f1=open('merge_26sample_filter_add_4method.txt','r')
f2=open('merge_26sample_filter_add_4method_count_no_merge_proportion.txt','w')
f2.write("Individual\tbase_change\tnumber\tproportion\n")

lines=f1.readlines()
f1.close()
types=[]
samples=[]
sampletype=[]
dict1={}
dict1_sum={}
for i in lines:
	tmp=i.strip().split('\t')
	types.append(tmp[1])
	samples.append(tmp[0])
	dict1.setdefault(tmp[0],[]).append(int(tmp[2]))
types=list(set(types))
samples=list(set(samples))
for i in samples:
	dict1_sum[i]=sum(dict1[i])
	for j in types:
		sampletype.append(i+'_'+j)
sampletype2=[]
for i in lines:
	tmp=i.strip().split('\t')
	pro=float(tmp[-1])/dict1_sum[tmp[0]]
	sampletype2.append(tmp[0]+'_'+tmp[1])
	f2.write(tmp[0]+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+str(round(pro,4))+'\n')
for i in sampletype:
	if i not in sampletype2:
		(a,b)=i.split('_')
		f2.write(a+'\t'+b+'\t0\t0'+'\n')
f2.close()

