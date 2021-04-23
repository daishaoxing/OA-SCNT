# -*- coding: utf-8 -*-
####usage:python step11_replace_genome.py 08431GFP
#初始传递命令

import os,sys,re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
path='/home/devdata/nyy/nyy_GFP/replace_genome'
genome='/home/devdata/Genome/fasta_gtf/Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa'
sample=sys.argv[1]
vcf ='../filter_vcf/'+sample+'.lofreq2_UMv1.vcf' ###need change
output_fasta = os.path.join(path,sample+'_Mmul_snp.fa')

##读入参考基因组序列
genome_seq_snp = {}
seq_id_all = []
fasta_id = []
 
for seq_record in SeqIO.parse(genome, "fasta"):
    seqname=seq_record.description
    #print (seq_record.description)
    m=re.search(r'(.+) dna\:',seqname)
    if len(m.groups()[0])<4:
        genome_seq_snp[m.groups()[0]]=list(seq_record.seq)
print (list(genome_seq_snp.keys()))

##读入 SNP 检测文件，替换参考基因组序列
with open(vcf, encoding='utf-8', errors='ignore') as snp_vcf:
	for line in snp_vcf:
		if line.strip()[0:2] != '##':
			tmp = line.strip().split('\t')
			if tmp[6]=='PASS' and len(tmp[3])==1 and len(tmp[4])==1:
				AF=tmp[7].split(';')[1][3:]
				if float(AF)>0.5:
					if tmp[0] in genome_seq_snp.keys():
						genome_seq_snp[tmp[0]][int(tmp[1]) - 1] = tmp[4]
snp_vcf.close()
 
#输出fasta
output_genome = open(output_fasta, 'w')
for i in genome_seq_snp.keys():
	output_genome.write('>'+i+'\n'+''.join(genome_seq_snp[i])+'\n')
output_genome.close()