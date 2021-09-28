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
f3a=open('all_SCNT_RNA_edit_siteseq.txt','w')
f3b=open('all_ABE_RNA_edit_siteseq.txt','w')
dict_site={}
for i in f1.readlines():
	tmp=i.strip().split('\t')
	(chrpos,strand)=tmp[1].split('|')
	(chrm,pos)=chrpos.split(':')
	pos1=int(pos)-1-1
	pos2=int(pos)+1
	seqname='>'+chrm
	subseq=dict1[chrm][pos1:pos2]
	if strand=='-':
		subseq=complement(subseq)
	if len(subseq)==3 and subseq[1]=='A':
		dict_site.setdefault(tmp[0],[]).append(chrm+'_'+pos+'\t'+subseq)
	else:
		print (seqname+'\n'+subseq)

for i in dict_site.keys():
	fname=i+'.RNA_edit_siteseq.txt'
	f2=open(fname,'w')
	f2.write('\n'.join(dict_site[i]))
	f2.close()
	if i in ["NC1","NC2","NC3","NC4","NC5"]:
		f3a.write('\n'.join(dict_site[i])+'\n')
	if i in ["PC1","PC2","PC3","PC4","PC5"]:
		f3b.write('\n'.join(dict_site[i])+'\n')

f1.close()
f3a.close()
f3b.close()