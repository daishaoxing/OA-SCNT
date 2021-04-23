import os

list1=['07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR']
# list1=['cloning-W','cloning-PF','cloning-JR']
###1	220	.	G	A	160	PASS	BRF=0.56;FR=0.5000;HP=2;HapScore=2;MGOF=53;MMLQ=37;MQ=60.0;NF=3;NR=3;PP=160;QD=32.4618687579;SC=TAAGACGTCCGTACCAGCGTG;SbPval=0.52;Source=Platypus;TC=7;TCF=3;TCR=4;TR=6;WE=228;WS=210
###	GT:GL:GOF:GQ:NR:NV	1/0:-20.09,0.0,-5.29:53:53:7:6

for i in list1:
	fname=i+'.platypus_UMv1.vcf'
	print (i)
	f1=open(fname,'r',errors='ignore')
	f2=open("../filter_vcfnew/"+i+'_platypus_SNP.txt','w')
	headinfo="POS\tREF\tALT\tGT\tAF\tAF2\tDP\tinfo\n"
	f2.write(headinfo)
	for i in f1.readlines():
		if i[0] !='#':
			tmp=i.strip().split('\t')
			if len(tmp[0])<3 and tmp[6]=='PASS' and len(tmp[4])==1 and len(tmp[3])==1:
				info=tmp[-1]              
				tmpinfo=info.split(':') ###INFO
				pos=tmp[0]+':'+tmp[1]
				ref=tmp[3]
				alt=tmp[4]
				GQ=float(tmpinfo[-3])
				dp=float(tmpinfo[-2])
				dp1=float(tmpinfo[-1])
				af=float(dp1/dp)
				af1=1-af
				if alt !='.' and af >0.10 and dp>=20 and GQ>30:
					if af1>0.10:
						gt=ref+alt
						f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[3][4:]])+'\n')
					else:
						gt=alt+alt
						f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[3][4:]])+'\n')
	f1.close()
	f2.close()
