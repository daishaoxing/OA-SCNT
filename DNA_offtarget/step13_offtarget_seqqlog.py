import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
dict1={}

for seq_record in SeqIO.parse("/home/devdata/nyy/nyy_GFP/replace_genome/08431GFP_Mmul_snp.fa", "fasta"):
	seqname=seq_record.description
#	print (seq_record.description)
	dict1[seqname]=str(seq_record.seq)
# print (list(dict1.keys()))

def complement(strA):
	dna_seq = Seq(strA)
	return str(dna_seq.reverse_complement())

f1=open('merge_15sample_filter_add_4method_overlapwithsamplename.txt','r')
f2=open('all_edit_site_seqlog.txt','w')
for i in f1.readlines():
	tmp=i.strip().split('\t')
	if tmp[1]=='A>G' or tmp[1]=='T>C':
		chr=tmp[0].split(':')[0]
		pos=tmp[0].split(':')[1]
		pos1=int(pos)-1-1
		pos2=int(pos)+1
		subseq=dict1[chr][pos1:pos2]
		if	tmp[1]=='A>G':
			if subseq[1]=='A':
				f2.write(chr+'_'+pos+'\t'+subseq+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+tmp[3]+'\n')
			else:
				print (i+subseq)
		if	tmp[1]=='T>C':
			subseqs=dict1[chr][pos1:pos2]
			subseq=complement(subseqs)
			print(subseq)
			if subseq[1]=='A':
				f2.write(chr+'_'+pos+'\t'+subseq+'\t'+tmp[1]+'\t'+tmp[2]+'\t'+tmp[3]+'\n')
			else:
				print (i+subseq)
f2.close()
f1.close()
