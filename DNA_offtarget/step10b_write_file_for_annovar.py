import os
f1=open('gatk_all_A2Gv2.txt','r')
f2=open('gatk_all_A2Gv2_snpfile_for_annovar.txt','w')
sample_dict={}
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')[2].split('|')
	(chr,pos)=tmp[0].split(':')
	(basea,baseb)=tmp[1].split('>')
	f2.write('\t'.join([chr,pos,pos,basea,baseb,'comments:',i.strip().split('\t')[2]])+'\n')
f1.close()
f2.close()
cmdstr='perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out '+'gatk_all_A2G_annovar_out -build rheMac8 '+'gatk_all_A2Gv2_snpfile_for_annovar.txt'+' /home/dsx/bin/annovar/rheMac8db -neargene 10000'
print(cmdstr)
os.system(cmdstr)
