import os,re
fout=open('all_pos_unique_sort.bed','w')
pos_list=[]
for i in os.listdir('/home/devdata/nyy/nyy_GFP/GFP/vcf_2021/'):
	if i.endswith('_SNP.txt'):
		fname=i
		f1=open(fname,'r')
		print(i)
		for j in f1.readlines()[1:]:
			tmp=j.strip().split('\t')
			pos_list.append(tmp[0])
		f1.close()
pos_list1=list(set(pos_list))
pos_sort=sorted(pos_list1)
print(pos_sort[0])
print(pos_sort[-1])
for i in pos_sort:
	tmp=i.strip().split(':')
	if len(tmp)==2:
		fout.write('\t'.join([tmp[0],tmp[1],tmp[1]])+'\n')
	else:
		print(i)
fout.close()


####222.197.214.27: /data/nyy/offtarget/add_fq/depth
####NC6.depth20.txt

#awk -F"\t" 'NR>1{print $1}' *_SNP.txt >all_pos.txt
# sort -u all_pos.txt >all_pos_unique_sort.txt
# awk -F":" '{print $1"\t"$2"\t"$2}' all_pos_unique_sort.txt >all_pos_unique.bed

#### awk '!a[$0]++' all_pos.txt >all_pos_unique.txt ####slow
######https://blog.csdn.net/feng973/article/details/73849586
#### sort all_pos_unique.txt >all_pos_unique_sort.txt
#### awk -F":" '{print $1"\t"$2"\t"$2}' all_pos_unique_sort.txt >all_pos_unique.bed
#### sambamba depth base --min-coverage 40 --regions=all_pos_unique.bed bamfile
