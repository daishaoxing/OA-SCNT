import os

list1=['SCNT-NC1', 'SCNT-NC2', 'SCNT-NC3', 'SCNT-NC4', 'SCNT-NC5','SCNT-PC1','SCNT-PC2','SCNT-PC3','SCNT-PC4','SCNT-PC5']

for i in list1:
	fname=i+'_gatk_filteredv2.vcf'
	print (i)
	f1=open(fname,'r',errors='ignore')
	f2=open(i+'_GATK_SNP05.txt','w')
	headinfo="POS\tREF\tALT\tGT\tAF\tAF2\tDP\n"
	f2.write(headinfo)
	for i in f1.readlines():
		if i[0] !='#':
			tmp=i.strip().split('\t') ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SCNT-PC1
			####GT:AD:DP:GQ:PL  1/1:0,8:8:24:255,24,0
			info=tmp[-1]
			tmpinfo=info.split(':') ###INFO
			# print (tmpinfo)
			if len(tmp[0])<3 and tmp[6]=='PASS' and len(tmpinfo)==5 and len(tmp[4])==1 and len(tmp[3])==1:  
				pos=tmp[0]+':'+tmp[1]
				ref=tmp[3]
				alt=tmp[4]
				dp=float(tmpinfo[2])
				if dp>0:
					(ad1,ad2)=tmpinfo[1].split(',')
					af=float(ad2)/dp
					af1=float(ad1)/dp
					if af >0.05:
						if af1>0.05:
							gt=ref+alt
							f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp)])+'\n')
						else:
							gt=alt+alt
							f2.write('\t'.join([pos,ref,alt,gt,str(af),str(af1),str(dp)])+'\n')
				else:
					print (tmpinfo)
	f1.close()
	f2.close()
