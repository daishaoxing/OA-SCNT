import os
f1=open('lofreq2_strelka2_overlap_filter_more_v1.txt','r')
f1a=open('genome_filter_snps_for_annovar8.txt','w')
f2=open('GATK_10sample_offtarget.txt','r')
f2a=open('RNAseq_filter_snps_for_annovar8.txt','w')
f3=open('genome_RNAseq_offtarget.txt','w')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')
	(chr,pos)=tmp[0].split(':')
	f1a.write('\t'.join([chr,pos,pos,tmp[2][1],tmp[6][1],'comments:none'])+'\n')
	f3.write(tmp[0]+'\tgenome\n')
for i in f2.readlines()[1:]:
	tmp=i.strip().split('\t')
	(chr,pos)=tmp[0].split(':')
	f2a.write('\t'.join([chr,pos,pos,tmp[2][1],tmp[6][1],'comments:none'])+'\n')
	f3.write(tmp[0]+'\ttranscriptome\n')

os.system("perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out genome_filter_more_snps -build rheMac8 genome_filter_more_snps_for_annovar8.txt /home/dsx/bin/annovar/rheMac8db -neargene 10000")
          #perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out filter_more_snps -build rheMac8 filter_more_snps_for_annovar8.txt /home/dsx/bin/annovar/rheMac8db -neargene 10000
os.system("perl /home/dsx/bin/annovar/annotate_variation.pl -geneanno -dbtype refGene -out RNAseq_filter_more_snpsv5 -build rheMac8 RNAseq_filter_v5_snps_for_annovar8.txt /home/dsx/bin/annovar/rheMac8db -neargene 10000")
