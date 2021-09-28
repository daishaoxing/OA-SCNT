import os

list1=['07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR']
# list1=['cloning-W','cloning-PF','cloning-JR']
for i in list1:
	fname=i+'.lofreq2_UMv1.vcf'
	print (i)
	f1=open(fname,'r',errors='ignore')
	f2=open("../filter_vcfnew/"+i+'_lofreq2_SNP.txt','w')
	headinfo="POS\tREF\tALT\tGT\tAF\tAF2\tDP\tinfo\n"
	f2.write(headinfo)
	for i in f1.readlines():
		if i[0] !='#':
			tmp=i.strip().split('\t') ##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
			####1       220     .       G       A       1690    PASS    AF=0.428571;DP=133;DP4=42,31,17,40;SB=26
			info=tmp[7]              
			tmpinfo=info.split(';') ###INFO
			if len(tmp[0])<3 and tmp[6]=='PASS' and len(tmp[4])==1 and len(tmp[3])==1:  
				pos=tmp[0]+':'+tmp[1]
				ref=tmp[3]
				alt=tmp[4]
				dp=float(tmpinfo[1][3:])
				af=float(tmpinfo[0][3:])
				dp4=tmpinfo[2][4:].split(',')
				af1=(float(dp4[0])+float(dp4[1]))/dp
				if alt !='.' and af >0.10 and dp>=20:
					if af1>0.10:
						gt=ref+alt
						f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[3][4:]])+'\n')
					else:
						gt=alt+alt
						f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp),tmpinfo[3][4:]])+'\n')
	f1.close()
	f2.close()
