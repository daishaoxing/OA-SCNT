import os
region=['exonic','intronic','intronic','UTR3','UTR5','ncRNA_exonic','ncRNA_intronic','UTR5;UTR3','UTR3;UTR5']
AG=['A>G','T>C']
fout=open('DNA_A2Ggene.txt','w')
fin=open('/home/devdata/nyy/nyy_GFP/GFP/filter_vcfnew/annovar_out.variant_function','r')
genedict={}
for j in fin.readlines():
	tmp=j.strip().split('\t')
	if tmp[0] in region:
		gene=tmp[1].split('.')[0]
		sample=tmp[-1].split('|')[0]
		change=tmp[-1].split('|')[1]
		if change in AG:
			genedict.setdefault(sample,[]).append(gene)
	else:
		print (tmp[0])
for j in genedict:
	genelist=list(set(genedict[j]))
	for i in genelist:
		fout.write(i+'\t'+j+'\n')
fout.close()