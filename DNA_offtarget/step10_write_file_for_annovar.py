import os
f1=open('merge_26sample_filter_add_4method_overlapwithsamplenamewithNC.txt','r')
f2=open('snpfile_for_annovar.txt','w')
sample_dict={}
for i in f1.readlines():
	tmp=i.strip().split('\t')
	(chr,pos)=tmp[0].split(':')
	(basea,baseb)=tmp[1].split('>')
	f2.write('\t'.join([chr,pos,pos,basea,baseb,'comments:',tmp[-1]+'|'+tmp[1]])+'\n')
f1.close()
f2.close()
cmdstr='perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out '+'annovar_out -build rheMac8 '+'snpfile_for_annovar.txt'+' /home/dsx/bin/annovar/rheMac8db -neargene 10000'
print(cmdstr)
os.system(cmdstr)
