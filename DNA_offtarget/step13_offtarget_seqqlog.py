import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
dict1={}

for seq_record in SeqIO.parse("/nas01/nyy/vcffrom27/replace_genome/08431GFP_Mmul_snp.fa", "fasta"):
	seqname=seq_record.description
#	print (seq_record.description)
	dict1[seqname]=str(seq_record.seq)
# print (list(dict1.keys()))

def complement(strA):
	dna_seq = Seq(strA)
	return str(dna_seq.reverse_complement())

f1=open('merge_26sample_filter_add_4method_overlapwithsamplenamewithNC.txt','r')
f2=open('all_edit_site_seqlog_AG.txt','w')
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
				f2.write(chr+'_'+pos+'\t'+subseq+'\t'+tmp[1]+'\t'+tmp[2]+'\n')
			else:
				print (i+subseq)
		if	tmp[1]=='T>C':
			subseqs=dict1[chr][pos1:pos2]
			subseq=complement(subseqs)
			if subseq[1]=='A':
				f2.write(chr+'_'+pos+'\t'+subseq+'\t'+tmp[1]+'\t'+tmp[2]+'\n')
			else:
				print (i+subseq)
f2.close()
f1.close()


f1=open('merge_26sample_filter_add_4method_overlapwithsamplenamewithNC.txt','r')
f2=open('all_edit_site_seqlog_CT.txt','w')
for i in f1.readlines():
	tmp=i.strip().split('\t')
	if tmp[1]=='G>A' or tmp[1]=='C>T':
		chr=tmp[0].split(':')[0]
		pos=tmp[0].split(':')[1]
		pos1=int(pos)-1-1
		pos2=int(pos)+1
		subseq=dict1[chr][pos1:pos2]
		if	tmp[1]=='C>T':
			if subseq[1]=='C':
				f2.write(chr+'_'+pos+'\t'+subseq+'\t'+tmp[1]+'\t'+tmp[2]+'\n')
			else:
				print (i+subseq)
		if	tmp[1]=='G>A':
			subseqs=dict1[chr][pos1:pos2]
			subseq=complement(subseqs)
			if subseq[1]=='C':
				f2.write(chr+'_'+pos+'\t'+subseq+'\t'+tmp[1]+'\t'+tmp[2]+'\n')
			else:
				print (i+subseq)
f2.close()
f1.close()

####Cytosine and adenine base editors (CBEs and ABEs), 
####which are fusions of a nickase-type Cas9 (nCas9) protein with a deaminase domain, 
####can catalyze the conversion of C→T (C>T) and A>G in the target site of a single guide RNA (sgRNA), respectively