fout=open('merge_26sample_change_raw.txt','w')

fname1='merge_26sample_change_raw_gatk.txt'
fname2='merge_26sample_change_raw_platypus.txt'
fname3='merge_26sample_change_raw_lofreq2.txt'
fname4='merge_26sample_change_raw_strelka2.txt'
gatk={}
platypus={}
lofreq={}
strelka={}
f1=open(fname1,'r',errors='ignore')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	gatk.setdefault(tmp[2],[]).append(tmp[0])
f1.close()
f1=open(fname2,'r',errors='ignore')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	platypus.setdefault(tmp[2],[]).append(tmp[0])
f1.close()
f1=open(fname3,'r',errors='ignore')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	lofreq.setdefault(tmp[2],[]).append(tmp[0])
f1.close()
f1=open(fname4,'r',errors='ignore')
for j in f1.readlines()[1:]:
	tmp=j.strip().split('\t')
	strelka.setdefault(tmp[2],[]).append(tmp[0])
f1.close()
for i in gatk.keys():
	overlap=list(set(gatk[i]).intersection(strelka[i],platypus[i],lofreq[i]))
	print(len(overlap))
	for j in overlap:
		fout.write(j+'\t'+i+'\n')
fout.close()