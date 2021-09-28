import os
f2=open('merge_26lofreq2_indel.txt','w')
headinfo="POS\ttype\tDP\tsample\n"
f2.write(headinfo)
for j in os.listdir('/data1/nyy/offtarget/v4indel'):
	if j.endswith('.lofreq2_UMv1.vcf'):
		samplename=j[0:-17]
		print (samplename)
		f1=open(j,'r',errors='ignore')
		for i in f1.readlines():
			if i[0] !='#':
				tmp=i.strip().split('\t') ##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
				####INFO		DP=20;AF=0.500000;SB=10;DP4=0,10,4,6
				####INFO		DP=24;AF=1.000000;SB=0;DP4=0,0,6,18;INDEL;HRUN=1
				info=tmp[-1]
				tmpinfo=info.split(';') ###INFO
				alt=tmp[4].split(',')
				if len(tmp[0])<3 and tmp[6]=='PASS' and (len(tmp[3])>1 or len(alt[0])>1) and len(alt)==1:
					pos=tmp[0]+':'+tmp[1]
					dp=float(tmpinfo[1][3:])
					af=float(tmpinfo[0][3:])
					if af >0.10 and dp>=40:
						change='Insertion'
						if len(tmp[3])>len(alt[0]):
							change='Deletion'
						f2.write('\t'.join([pos,change,str(dp),samplename])+'\n')
		f1.close()
f2.close()
