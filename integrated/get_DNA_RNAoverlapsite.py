import os
f1=open('GATK_10sample_edit_site.txt','r')
f2=open('lofreq2_strelka2_overlap_filter_more_edit_siteAG_CT.txt','r')
f3=open('genome_RNAseq_overlapAG_CT.txt','w')
genome=[]
RNAseq=[]
for i in f1.readlines():
	tmp=i.strip().split('\t')
	if tmp[0][0:2]=='PC':
		RNAseq.extend(tmp[2:])
		
for i in f2.readlines():
	tmp=i.strip().split('\t')
	if tmp[0][0:2]=='PC':
		genome.extend(tmp[2:])
genome=list(set(genome))
RNAseq=list(set(RNAseq))
f3.write("DNA\n"+"\n".join(genome))
f3.write("\nRNA\n"+"\n".join(RNAseq))