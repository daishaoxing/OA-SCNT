import os,re

list1=['SCNT-NC1', 'SCNT-NC2', 'SCNT-NC3', 'SCNT-NC4', 'SCNT-NC5','SCNT-PC1','SCNT-PC2','SCNT-PC3','SCNT-PC4','SCNT-PC5']

for j in list1:
	fname=j+'_gatk_filtered.vcf'
	print (j)
	f1=open(fname,'r',errors='ignore')
	f2=open(j+'_gatk_filteredv2.vcf','w')
	for i in f1.readlines():
		if i[0] !='#':
			tmp=i.strip().split('\t') 
			##CHROM	 POS	   ID	   REF	   ALT	   QUAL	   FILTER  INFO	   FORMAT  SCNT-PC1
			####GT:AD:DP:GQ:PL	1/1:0,8:8:24:255,24,0
			info=tmp[-1]
			tmpinfo=info.split(':') ###INFO
			MQ=0
			m=re.search(r'MQ=(.+?);',tmp[7])
			if m:
				MQ=float(m.groups()[0])
			if len(tmp[0])<3 and tmp[6]=='PASS' and len(tmpinfo)==5 and len(tmp[4])==1 and len(tmp[3])==1 and float(tmp[5])>25 and MQ>20:
				dp=float(tmpinfo[2])
				if dp>20:
					f2.write(i)
	f1.close()
	f2.close()
