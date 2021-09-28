import os
f1=open('merge_GATK_10sample_Jitter.txt','r')
dict_site={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t') ###PC1     25.0    PC      19      19:50150877
	(chr,pos)=tmp[4].split(':')
	dict_site.setdefault(tmp[0],[]).append([chr,pos,tmp[1],tmp[-1]])
f1.close()
# print(dict_site)
for i in ['PC1','PC2','PC3','PC4','PC5']:
	fname=i+'_snps_for_annovar8.txt'
	f2=open(fname,'w')
	for j in dict_site[i]:
		(chr,pos,freq,ref_alt)=j
		f2.write('\t'.join([chr,pos,pos,'A','G','comments:|'+freq+'|'+ref_alt])+'\n')
	f2.close()
	cmdstr='perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out '+i+'_RNA_filter_snps -build rheMac8 '+fname+' /home/dsx/bin/annovar/rheMac8db -neargene 10000'
	os.system(cmdstr)
