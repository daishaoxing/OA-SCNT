import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
dict1={}

for seq_record in SeqIO.parse("/home/devdata/nyy/nyy_GFP/replace_genome/08431GFP_Mmul_snp.fa", "fasta"):
	seqname=seq_record.description
	dict1[seqname]=str(seq_record.seq)
print (list(dict1.keys()))

f1=open('selected_A2G_more.txt','r')
f2=open('selected_A2G_offtarget_more_seq_2000.fa','w')
for i in f1.readlines()[1:]:
	tmp=i.strip().split('\t')[4].split('|')[0].split(':')
	chr=tmp[0]
	pos=tmp[1]
	pos1=int(tmp[1])-1000-1
	pos2=int(tmp[1])+1000
	seqname='>'+chr+':'+str(pos1)+'_'+pos+'_'+str(pos2)
	subseq=dict1[chr][pos1:pos2]
	f2.write(seqname+'\n'+subseq+'\n')
f1.close()
f2.close()
