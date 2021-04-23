import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
dict1={}

for seq_record in SeqIO.parse("/home/devdata/nyy/nyy_GFP/replace_genome/08431GFP_Mmul_snp.fa", "fasta"):
	seqname=seq_record.description
	print (seq_record.description)
	dict1[seqname]=str(seq_record.seq)
print (list(dict1.keys()))

def complement(strA):
	dna_seq = Seq(strA)
	return str(dna_seq.reverse_complement())

f1=open('GATK_10sample_edit_site.txt','r')

allsite=[]
f3=open('all_RNA_edit_siteseq.txt','w')
for i in f1.readlines():
	tmp=i.strip().split('\t')
	fname=tmp[0]+'.RNA_edit_siteseq.txt'
	f2=open(fname,'w')
	for j in tmp[2:]:
		info=j.split('|')
		strand=info[1]
		chr=info[0].split(':')[0]
		pos=info[0].split(':')[1]
		pos1=int(pos)-1-1
		pos2=int(pos)+1
		seqname='>'+info[0]
		subseq=dict1[chr][pos1:pos2]
		if strand=='-':
			subseqs=dict1[chr][pos1:pos2]
			subseq=complement(subseqs)
			if subseq[1]=='A':
				f2.write(chr+'_'+pos+'\t'+subseq+'\n')
				allsite.append(chr+'_'+pos+'\t'+subseq+'\n')
			else:
				print (seqname+'\n'+subseq)
	f2.close()
f1.close()

for i in set(allsite):
	f3.write(i)