import os
# list1=['07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR']
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
##1       220     .       G       A       127     PASS    MQ=60;SNVHPOL=3 GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL    0/1:151:13:20:0:10,10:6,6:4,4:-17.3:PASS:162,0,149
f2=open('merge_26gatk_indel.txt','w')
headinfo="POS\ttype\tDP\tsample\n"
f2.write(headinfo)
for j in os.listdir('/data1/nyy/offtarget/v4indel'):
	if j.endswith('.gatk_indel_UMv1.vcf'):
		samplename=j[0:-20]
		print (samplename)
		f1=open(j,'r',errors='ignore')
		for i in f1.readlines():
			if i[0] !='#':
				tmp=i.strip().split('\t')
				alt=tmp[4].split(',')
				if len(tmp[0])<3 and tmp[6]=='PASS' and (len(tmp[3]) >1 or len(alt[0])>1)  and len(alt)==1:
					info=tmp[-1]
					pos=tmp[0]+':'+tmp[1]
					tmpinfo=info.split(':')
					dp=float(tmpinfo[2])
					GQ=float(tmpinfo[3])
					if len(alt)==1 and GQ>30 and dp>=40:
						AD=tmpinfo[1].split(',')
						af=float(AD[1])/dp
						if af >0.10:
							change='Insertion'
							if len(tmp[3])>len(alt[0]):
								change='Deletion'
							f2.write('\t'.join([pos,change,str(dp),samplename])+'\n')
		f1.close()
f2.close()
