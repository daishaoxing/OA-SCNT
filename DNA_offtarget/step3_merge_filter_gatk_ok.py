#############background SNP
population_pos=[]
f1=open('/nas01/nyy/vcffrom27/population_all_SNP.txt','r')
for i in f1.readlines():
	tmp=i.strip().split('\t')
	change=tmp[0]
	population_pos.append(change)
f1.close()

f1=open('/nas01/nyy/vcffrom27/population_all_SNP_p2.txt','r')
for i in f1.readlines():
	tmp=i.strip().split('\t')
	change=tmp[0]
	population_pos.append(change)
f1.close()

f1=open('/nas01/nyy/vcffrom27/population_all_SNP_p3.txt','r')
for i in f1.readlines():
	tmp=i.strip().split('\t')
	change=tmp[0]
	population_pos.append(change)
f1.close()

#############start filter background SNP
f1=open('merge_26sample_refgenome_change_gatk.txt','r')
f2=open('merge_26sample_adddonor_change_gatk.txt','r')
f3=open('merge_26sample_refgenome_adddonor_count_gatk.txt','w')
f4=open('merge_26sample_refgenome_adddonor_change_gatk.txt','w')
f3.write("Individual\tbase_change\tnumber\n")
f4.write("base_change\tIndividual\tpos\n")

dict1={}
dict2={}
all_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t') ####4:30262373      A>T     AT      NC4
	change=tmp[0]+'|'+tmp[1]
	dict1.setdefault(tmp[3],[]).append(change)
	all_dict.setdefault(tmp[3],[]).append(change)
f1.close()
for i in f2.readlines()[1:]:
	tmp=i.strip().split('\t') ####4:30262373      A>T     AT      NC4
	change=tmp[0]+'|'+tmp[1]
	dict2.setdefault(tmp[3],[]).append(change)
	all_dict.setdefault(tmp[3],[]).append(change)
f2.close()

final_dict={}
for i in all_dict.keys():
	tmplist=all_dict[i]
	print(i)
	filterlist=list(set(tmplist)-set(population_pos))
	for j in filterlist:
		(a,b)=j.split('|')
		final_dict.setdefault(i,[]).append(b)
		f4.write('\t'.join([a,b,i])+'\n')

for i in final_dict.keys():
	tmplist=final_dict[i]
	for j in list(set(tmplist)):
		f3.write(i+'\t'+j+'\t'+str(tmplist.count(j))+'\n')
f3.close()
f4.close()

####################only for test
f3a=open('merge_26sample_refgenome_count_gatk.txt','w')
f3b=open('merge_26sample_adddonor_count_gatk.txt','w')
f3a.write("Individual\tbase_change\tnumber\n")
f3b.write("Individual\tbase_change\tnumber\n")

final_dict1={}
for i in dict1.keys():
	tmplist=dict1[i]
	print(i)
	filterlist=list(set(tmplist)-set(population_pos))
	for j in filterlist:
		(a,b)=j.split('|')
		final_dict1.setdefault(i,[]).append(b)
for i in final_dict1.keys():
	tmplist=final_dict1[i]
	for j in list(set(tmplist)):
		f3a.write(i+'\t'+j+'\t'+str(tmplist.count(j))+'\n')

final_dict2={}
for i in dict2.keys():
	tmplist=dict2[i]
	print(i)
	filterlist=list(set(tmplist)-set(population_pos))
	for j in filterlist:
		(a,b)=j.split('|')
		final_dict2.setdefault(i,[]).append(b)
for i in final_dict2.keys():
	tmplist=final_dict2[i]
	for j in list(set(tmplist)):
		f3b.write(i+'\t'+j+'\t'+str(tmplist.count(j))+'\n')

f3a.close()
f3b.close()
