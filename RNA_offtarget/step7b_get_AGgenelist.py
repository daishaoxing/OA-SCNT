import os
list1=['PC1','PC2','PC3','PC4','PC5'] ###PC5_RNA_filter_snps.variant_function
region=['exonic','intronic','intronic','UTR3','UTR5','ncRNA_exonic','ncRNA_intronic','UTR5;UTR3','UTR3;UTR5','splicing', 'ncRNA_splicing', 'exonic;splicing']
fout1=open('A2Ggene.txt','w')
fout2=open('A2Ggene_freq.txt','w')
region_tmp=[]
for i in list1:
	fin=open(i+'_RNA_filter_snps.variant_function','r')
	genelist=[]
	for j in fin.readlines():
		tmp=j.strip().split('\t')
		if tmp[0] in region:
			comments=tmp[-1].split('|') ##comments:|20.369999999999997|11:848308|A>G
			pos_alt=comments[2]+'|'+comments[3]
			freq=comments[1]
			for k in tmp[1].split(';'):
				gene=k.split('.')[0]
				genelist.append(gene)
				fout2.write('\t'.join([i,gene,freq,tmp[2],tmp[3],pos_alt])+'\n')
		else:
			region_tmp.append(tmp[0])
	genelist=list(set(genelist))
	for j in genelist:
		fout1.write(j+'\t'+i+'\n')
print(set(region_tmp))
fout1.close()
fout2.close()