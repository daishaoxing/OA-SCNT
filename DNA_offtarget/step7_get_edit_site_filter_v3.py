import numpy as np
import _pickle as cPickle
f1=open('4method_overlap_site_freq.txt','r')
f2=open('gatk_all_A2Gv2.bed','w')
f3=open('gatk_all_A2Gv2.txt','w')
lines=f1.readlines()
f3.write('number\tsamples\t'+lines[0].strip()+'\tdiffrq1\tdiffrq2\n')
# AG=['A>G','T>C','C>T','G>A']
AG=['A>G','T>C']
PCsample=['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','cloningW','cloningPF','cloningJR']
depth40final = cPickle.load(open("depth40_pos.pkl","rb"))
def compare(a,b):
	diffrq=0
	if a !='NA':
		(GT1,freq1)=a.split('|')
		gtlist=[]
		frqlist=[]
		for i in b:
			if i !='NA':
				(GT2,freq2)=i.split('|')
				gtlist.append(GT2)
				frqlist.append(float(freq2))
		# print(gtlist)
		if GT1 not in gtlist and (len(set(gtlist))==1):
			if float(freq1)>np.mean(frqlist):
				diffrq=float(freq1)-np.mean(frqlist)
			else:
				diffrq=float(freq1)
	return diffrq

for i in lines[1:]:
	tmp=i.strip().split('\t')
	site=tmp[0].split('|')
	if site[0] in depth40final:
		if site[1] in AG:
			control=tmp[2:6]
			NAnum=control.count('NA')
			print(NAnum)
			edit=tmp[6:]
			if NAnum<4:
				oksample=[]
				freq=[]
				for j in range(len(edit)):
					diffrq=compare(edit[j],control)
					if diffrq>0.2:#####0.2 for subset of bam ,0.4 for filter more sites 
						oksample.append(PCsample[j])
						freq.append(str(diffrq))
				if len(oksample)>=1:
					(chr,pos)=site[0].split(':')
					pos1=str(int(pos)-150)
					pos2=str(int(pos)+150)
					f2.write('\t'.join([str(chr),pos1,pos2])+'\n')
					f3.write(str(len(site)-3)+"\t"+'|'.join(oksample)+'\t'+i.strip()+'\t'+'\t'.join(freq)+'\n')
f2.close()
f3.close()


# position|change|GT|num	refbase	GFP	NC1	NC2	NC3	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	cloningW	cloningPF	cloningJR
# 10:1001420|T>C|CC|3	C	TT|1.0	TT|1.0	TT|1.0	TT|1.0	NA	NA	NA	NA	NA	NA	TT|1.0	NA	NA	TT|1.0	TT|0.98
# 10:10186236|T>G|TG|2	T	NA	NA	NA	NA	NA	NA	TG|0.58	NA	TG|0.27	NA	NA	NA	NA	NA	NA
# 10:10186261|G>T|GT|2	G	NA	NA	NA	NA	NA	NA	GT|0.52	NA	GT|0.29	NA	NA	NA	NA	NA	NA
# 10:1025098|C>T|TT|2	T	TC|0.39	TC|0.48	TC|0.41	TC|0.53	NA	NA	NA	NA	TC|0.69	NA	TC|0.37	NA	TC|0.49	TC|0.49	TC|0.49
# 10:1031642|A>G|GG|1	G	GA|0.53	GA|0.45	GA|0.43	GA|0.63	GA|0.69	NA	NA	GA|0.35	NA	NA	GA|0.7	NA	GA|0.46	GA|0.68	GA|0.49
# 10:1032522|A>G|GG|1	G	GA|0.46	GA|0.5	GA|0.63	GA|0.47	GA|0.5	NA	NA	GA|0.19	GA|0.64	GA|0.36	GA|0.56	NA	GA|0.58	GA|0.46	GA|0.5
# 10:1032539|T>C|CC|1	C	CT|0.57	CT|0.53	CT|0.39	CT|0.5	CT|0.46	NA	NA	CT|0.77	CT|0.38	NA	CT|0.4	NA	CT|0.42	CT|0.5	CT|0.51
# 10:1033848|G>A|AA|1	A	AG|0.54	AG|0.5	AG|0.7	AG|0.55	AG|0.55	AG|0.34	NA	AG|0.34	AG|0.56	AG|0.61	AG|0.63	NA	AG|0.5	AG|0.47	AG|0.64
# 10:1034264|C>G|GG|1	G	GC|0.51	GC|0.54	GC|0.5	GC|0.4	GC|0.55	GC|0.36	NA	GC|0.39	GC|0.54	GC|0.5	GC|0.39	NA	GC|0.54	GC|0.53	GC|0.35
