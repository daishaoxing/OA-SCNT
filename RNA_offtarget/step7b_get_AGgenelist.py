import os
list1=['PC1','PC2','PC3','PC4','PC5'] ###PC5_RNA_filter_snps.variant_function
region=['exonic','intronic','intronic','UTR3','UTR5','ncRNA_exonic','ncRNA_intronic','UTR5;UTR3','UTR3;UTR5']
fout=open('A2Ggene.txt','w')
for i in list1:
	fin=open(i+'_RNA_filter_snps.variant_function','r')
	genelist=[]
	for j in fin.readlines():
		tmp=j.strip().split('\t')
		if tmp[0] in region:
			for k in tmp[1].split(';'):
				genelist.append(k.split('.')[0])
		else:
			print (j)
	genelist=list(set(genelist))
	for j in genelist:
		fout.write(j+'\t'+i+'\n')
fout.close()