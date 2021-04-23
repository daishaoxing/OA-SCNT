f2=open('SNP_gene_region_count.txt','w')
f2.write('sample\ttype\tRegion\tcount\n') ##07027\tSNP	Exonic	55974
fname="annovar_out.variant_function" ##
f1=open(fname,'r',errors='ignore')
dict1={}
list1=['Exonic','UTR','Intronic','Intergenic','Downstream','Upstream']
# ['UTR3', 'UTR5', 'downstream', 'exonic', 'intergenic', 'intronic', 'ncRNA_exonic', 'ncRNA_intronic', 'splicing', 'upstream', 'upstream;downstream']
for j in f1.readlines():
	tmp=j.strip().split('\t')
	mykey='\t'.join(tmp[-1].split('|'))
	if tmp[0] in ['UTR3','UTR5','UTR5;UTR3','UTR3;UTR5']:
		dict1.setdefault(mykey,[]).append('UTR')
	elif tmp[0] in ['exonic','ncRNA_exonic','splicing','ncRNA_splicing','exonic;splicing','ncRNA_exonic;splicing']:
		dict1.setdefault(mykey,[]).append('Exonic')
	elif tmp[0] in ['intronic','ncRNA_intronic']:
		dict1.setdefault(mykey,[]).append('Intronic')
	elif tmp[0] in ['intergenic']:
		dict1.setdefault(mykey,[]).append('Intergenic')
	elif tmp[0] in ['downstream','downstream;upstream']:
		dict1.setdefault(mykey,[]).append('Downstream')
	elif tmp[0] in ['upstream','upstream;downstream']:
		dict1.setdefault(mykey,[]).append('Upstream')
	else:
		print (j)
for i in dict1.keys():
	for j in set(dict1[i]):
		f2.write(i+'\t'+j+'\t'+str(dict1[i].count(j))+'\n')
f1.close()

