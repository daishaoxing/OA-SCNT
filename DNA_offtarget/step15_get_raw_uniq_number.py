import os
fout=open('raw_uniq_number.txt','w')
samplelist=['07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR']
for i in samplelist:
	fname1=i+'_4method_raw_dp20_filter_uniq.txt'
	f1=open(fname1,'r')
	fout.write(i+'\t'+str(len(f1.readlines()))+'\n')
