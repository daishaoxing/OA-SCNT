import os
f1=open('GATK_10sample_edit_site.txt','r')
for i in f1.readlines()[0:5]:
	tmp=i.strip().split('\t')
	fname=tmp[0]+'_snps_for_annovar8.txt'
	f2=open(fname,'w')
	for j in tmp[2:]:
		(chr,pos)=j[0:-2].split(':')
		f2.write('\t'.join([chr,pos,pos,'A','G','comments:none'])+'\n')
	f2.close()
	cmdstr='perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out '+tmp[0]+'_RNA_filter_snps -build rheMac8 '+fname+' /home/dsx/bin/annovar/rheMac8db -neargene 10000'
	os.system(cmdstr)
f1.close()