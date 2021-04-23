import os
AG=['A>G','T>C']
fout=open('DNA_A2Gsite.txt','w')
fin=open('/home/devdata/nyy/nyy_GFP/GFP/filter_vcfnew/snpfile_for_annovar.txt','r')
sitedict={}
for j in fin.readlines():
	tmp=j.strip().split('\t')
	site=':'.join([tmp[0],tmp[1]])
	sample=tmp[-1].split('|')[0]
	change=tmp[-1].split('|')[1]
	sitechange=site+'|'+change
	if change in AG:
		sitedict.setdefault(sample,[]).append(sitechange)
for j in sitedict:
	sitelist=list(set(sitedict[j]))
	for i in sitelist:
		fout.write(i+'\t'+j+'\n')
fout.close()