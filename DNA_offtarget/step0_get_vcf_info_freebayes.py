import os

list1=['07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR']
##GT:DP:AD:RO:QR:AO:QA:GL 0/1:142:76,60:76:2698:60:2177:-155.109,0,-201.97
for i in list1:
	fname=i+'.freebayes_UMv1.vcf'
	print (i)
	f1=open(fname,'r',errors='ignore')
	f2=open("../filter_vcfnew/"+i+'_freebayes_SNP.txt','w')
	headinfo="POS\tREF\tALT\tGT\tAF\tAF2\tDP\tinfo\n"
	f2.write(headinfo)
	for i in f1.readlines():
		if i[0] !='#':
			tmp=i.strip().split('\t')
			if len(tmp[0])<3 and len(tmp[4])==1 and len(tmp[3])==1:
				info=tmp[-1]              
				tmpinfo=info.split(':')
				pos=tmp[0]+':'+tmp[1]
				ref=tmp[3]
				alt=tmp[4]
				dp=float(tmpinfo[1])
				dp2=tmpinfo[2].split(',')
				af=float(dp2[0])/dp
				af1=float(dp2[1])/dp
				if af >0.10 and dp>=20:
					if af1>0.10:
						gt=ref+alt
						f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[3][4:]])+'\n')
					else:
						gt=alt+alt
						f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[3][4:]])+'\n')
	f1.close()
	f2.close()
