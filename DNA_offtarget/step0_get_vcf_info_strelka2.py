import os
import os
list1=['07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR']
# list1=['cloning-W','cloning-PF','cloning-JR']
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
##1       220     .       G       A       127     PASS    MQ=60;SNVHPOL=3 GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL    0/1:151:13:20:0:10,10:6,6:4,4:-17.3:PASS:162,0,149
for i in list1:
	fname=i+'.strelka2_UMv1.vcf'
	print (i)
	f1=open(fname,'r',errors='ignore')
	f2=open("../filter_vcfnew/"+i+'_strelka2_SNP.txt','w')
	headinfo="POS\tREF\tALT\tGT\tAF\tAF2\tDP\tinfo\n"
	f2.write(headinfo)
	for i in f1.readlines():
		if i[0] !='#':
			tmp=i.strip().split('\t')
			if len(tmp[0])<3 and tmp[6]=='PASS' and len(tmp[4])==1 and len(tmp[3])==1:
				info=tmp[-1]
				pos=tmp[0]+':'+tmp[1]
				ref=tmp[3]
				alt=tmp[4]
				tmpinfo=info.split(':')
				dp=float(tmpinfo[3])
				GQ=float(tmpinfo[1])
				if alt !='.' and GQ>30 and dp>=20:
					AD=tmpinfo[5].split(',')
					af=float(AD[1])/dp
					af1=float(AD[0])/dp
					if af >0.10:
						if af1>0.10:
							gt=ref+alt
							f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[5]])+'\n')
						else:
							gt=alt+alt
							f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[5]])+'\n')
	f1.close()
	f2.close()
