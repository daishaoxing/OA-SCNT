f1=open('all_pos.txt','r')
f2=open('all_pos_unique.bed','w')
list1=f1.readlines()
list1=list(set(list1))
list2=sorted(list1)
for i in list2:
	tmp=i.strip().split(':')
	if len(tmp)==2:
		f2.write('\t'.join([tmp[0],tmp[1],tmp[1]])+'\n')
	else:
		print(i)