import os
f2=open('merge_26platypus_indel.txt','w')
headinfo="POS\ttype\tDP\tsample\n"
f2.write(headinfo)
###1	220	.	G	A	160	PASS	BRF=0.56;FR=0.5000;HP=2;HapScore=2;MGOF=53;MMLQ=37;MQ=60.0;NF=3;NR=3;PP=160;QD=32.4618687579;SC=TAAGACGTCCGTACCAGCGTG;SbPval=0.52;Source=Platypus;TC=7;TCF=3;TCR=4;TR=6;WE=228;WS=210
###	GT:GL:GOF:GQ:NR:NV	1/0:-20.09,0.0,-5.29:53:53:7:6

for j in os.listdir('/data1/nyy/offtarget/v4indel'):
	if j.endswith('.platypus_UMv1.vcf'):
		samplename=j[0:-18]
		print (samplename)
		f1=open(j,'r',errors='ignore')
		for i in f1.readlines():
			if i[0] !='#':
				tmp=i.strip().split('\t')
				alt=tmp[4].split(',')
				if len(tmp[0])<3 and tmp[6]=='PASS' and (len(alt[0])>1 or len(tmp[3])>1) and len(alt)==1:
					info=tmp[-1]
					tmpinfo=info.split(':') ###INFO
					pos=tmp[0]+':'+tmp[1]
					GQ=float(tmpinfo[-3])
					dp=float(tmpinfo[-2])
					dp1=float(tmpinfo[-1])
					if dp>=40 and GQ>30:
						af=float(dp1/dp)
						if af >0.10:
							change='Insertion'
							if len(tmp[3])>len(alt[0]):
								change='Deletion'
							f2.write('\t'.join([pos,change,str(dp),samplename])+'\n')
		f1.close()
f2.close()
