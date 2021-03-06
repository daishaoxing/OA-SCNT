import os
f1=open('gatk_all_A2Gv2.txt','r')
f2=open('runIGV_batch_gatk_all_A2Gv2.txt','w')
f2.write('new\n')
f2.write('genome D:/kmust/nyy/genome_ref/Macaca_mulatta_revised.fa\n')

f2.write('load D:/kmust/nyy/genome_ref/08431GFP_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/NC1_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/NC2_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/NC3_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC1_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC2_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC3_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC4_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC5_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC6_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC7_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/PC8_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/cloning-JR_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/cloning-PF_offtarget.bam\n')
f2.write('load D:/kmust/nyy/genome_ref/cloning-W_offtarget.bam\n')
f2.write('snapshotDirectory D:/kmust/nyy/genome_ref\n')
num=1
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	info=tmp[2].split('|')
	(chr,pos)=info[0].split(':')
	pos1=str(int(pos)-20)
	pos2=str(int(pos)+20)
	str1=('goto chr'+chr+':'+pos1+'-'+pos2+'\n')
	f2.write(str1)
	f2.write('sort base\n')
	f2.write('squish\n')
	f2.write('snapshot selected_position'+'_'.join(info[3:])+'_'+str(num)+'.png\n')
	num=num+1
