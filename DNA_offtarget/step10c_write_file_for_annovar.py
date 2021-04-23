import os
f1=open('PCR_site.txt','r')
f2=open('PCR_site_for_annovar.txt','w')
sample_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')[5].split('|')
	(chr,pos)=tmp[0].split(':')
	(basea,baseb)=tmp[1].split('>')
	f2.write('\t'.join([chr,pos,pos,basea,baseb,'comments:',i.strip().split('\t')[5]])+'\n')
f1.close()
f2.close()
cmdstr='perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out '+'PCR_site_annovar_out -build rheMac8 '+'PCR_site_for_annovar.txt'+' /home/dsx/bin/annovar/rheMac8db -neargene 10000'
print(cmdstr)
os.system(cmdstr)
