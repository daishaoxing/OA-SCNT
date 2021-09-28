f2=open('SNP_gene_region_count_all.txt','w')
f3=open('SNP_gene_region_count_AG.txt','w')
f4=open('SNP_gene_region_count_CT.txt','w')
f2.write('sample\tRegion\tnumber\n') ##07027\tSNP	Exonic	55974
fname="annovar_out.variant_function" ##
f1=open(fname,'r',errors='ignore')
lines=f1.readlines()
f1.close()

list1=['Exonic','UTR','Intronic','Intergenic','Downstream','Upstream']
# ['UTR3', 'UTR5', 'downstream', 'exonic', 'intergenic', 'intronic', 'ncRNA_exonic', 'ncRNA_intronic', 'splicing', 'upstream', 'upstream;downstream'
sample_list=[]
sample_region=[]
for j in lines:
	tmp=j.strip().split('\t')
	sample=tmp[-1].split('|')[0]
	sample_list.append(sample)
	if tmp[0] in ['UTR3','UTR5','UTR5;UTR3','UTR3;UTR5']:
		sample_region.append(sample+'\tUTR')
	elif tmp[0] in ['exonic','ncRNA_exonic','splicing','ncRNA_splicing','exonic;splicing','ncRNA_exonic;splicing']:
		sample_region.append(sample+'\tExonic')
	elif tmp[0] in ['intronic','ncRNA_intronic']:
		sample_region.append(sample+'\tIntronic')
	elif tmp[0] in ['intergenic']:
		sample_region.append(sample+'\tIntergenic')
	elif tmp[0] in ['downstream','downstream;upstream']:
		sample_region.append(sample+'\tDownstream')
	elif tmp[0] in ['upstream','upstream;downstream']:
		sample_region.append(sample+'\tUpstream')
	else:
		print (j)
sample_region2=[]
for i in set(sample_list):
	for j in list1:
		sample_region2.append(i+'\t'+j)

for i in sample_region2:
	if i in sample_region:
		f2.write(i+'\t'+str(sample_region.count(i))+'\n')
	else:
		f2.write(i+'\t0\n')
f2.close()



##########c('A>G','T>C'),]
f3.write('sample\tRegion\tnumber\n')
list1=['Exonic','UTR','Intronic','Intergenic','Downstream','Upstream']
sample_list=[]
sample_region=[]
for j in lines:
	tmp=j.strip().split('\t')
	sample=tmp[-1].split('|')[0]
	change=tmp[-1].split('|')[1]
	sample_list.append(sample)
	if change in ['A>G','T>C']:
		if tmp[0] in ['UTR3','UTR5','UTR5;UTR3','UTR3;UTR5']:
			sample_region.append(sample+'\tUTR')
		elif tmp[0] in ['exonic','ncRNA_exonic','splicing','ncRNA_splicing','exonic;splicing','ncRNA_exonic;splicing']:
			sample_region.append(sample+'\tExonic')
		elif tmp[0] in ['intronic','ncRNA_intronic']:
			sample_region.append(sample+'\tIntronic')
		elif tmp[0] in ['intergenic']:
			sample_region.append(sample+'\tIntergenic')
		elif tmp[0] in ['downstream','downstream;upstream']:
			sample_region.append(sample+'\tDownstream')
		elif tmp[0] in ['upstream','upstream;downstream']:
			sample_region.append(sample+'\tUpstream')
		else:
			print (j)
sample_region2=[]
for i in set(sample_list):
	for j in list1:
		sample_region2.append(i+'\t'+j)

for i in sample_region2:
	if i in sample_region:
		f3.write(i+'\t'+str(sample_region.count(i))+'\n')
	else:
		f3.write(i+'\t0\n')
f3.close()


#########################c('C>T','G>A')
f4.write('sample\tRegion\tnumber\n')
list1=['Exonic','UTR','Intronic','Intergenic','Downstream','Upstream']
sample_list=[]
sample_region=[]
for j in lines:
	tmp=j.strip().split('\t')
	sample=tmp[-1].split('|')[0]
	change=tmp[-1].split('|')[1]
	sample_list.append(sample)
	if change in ['C>T','G>A']:
		if tmp[0] in ['UTR3','UTR5','UTR5;UTR3','UTR3;UTR5']:
			sample_region.append(sample+'\tUTR')
		elif tmp[0] in ['exonic','ncRNA_exonic','splicing','ncRNA_splicing','exonic;splicing','ncRNA_exonic;splicing']:
			sample_region.append(sample+'\tExonic')
		elif tmp[0] in ['intronic','ncRNA_intronic']:
			sample_region.append(sample+'\tIntronic')
		elif tmp[0] in ['intergenic']:
			sample_region.append(sample+'\tIntergenic')
		elif tmp[0] in ['downstream','downstream;upstream']:
			sample_region.append(sample+'\tDownstream')
		elif tmp[0] in ['upstream','upstream;downstream']:
			sample_region.append(sample+'\tUpstream')
		else:
			print (j)
sample_region2=[]
for i in set(sample_list):
	for j in list1:
		sample_region2.append(i+'\t'+j)

for i in sample_region2:
	if i in sample_region:
		f4.write(i+'\t'+str(sample_region.count(i))+'\n')
	else:
		f4.write(i+'\t0\n')
f3.close()