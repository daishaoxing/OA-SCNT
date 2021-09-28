import os
import _pickle as cPickle

f1=open('merge_26platypus_indel.txt','r')
f2=open('merge_26platypus_filter_indel.txt','w')
depth40list = cPickle.load(open("/nas01/nyy/vcffrom27/depth40_vcfpos.pkl","rb"))
GFP=['08431neg','08431pos','08431GFP','08431WT']
SCNT=['090202body','090202head','NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9']
ABE=['cloning-JR','cloning-PF','cloning-W','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8']

Ins_dict={}
del_dict={}
for i in f1.readlines()[1:]:
    tmp=i.strip().split('\t')##1:9580  Insertion       25.0    07027
    if tmp[1]=='Insertion' and float(tmp[2])>=40 and tmp[3] !='07027':
        Ins_dict.setdefault(tmp[3],[]).append(tmp[0])
    if tmp[1]=='Deletion' and float(tmp[2])>=40 and tmp[3] !='07027':
        del_dict.setdefault(tmp[3],[]).append(tmp[0])
####Insertion
GFP_list=[]
for i in GFP:
    GFP_list=GFP_list+Ins_dict[i]
SCNT_list=[]
for i in SCNT:
    SCNT_list=SCNT_list+Ins_dict[i]
ABE_list=[]
for i in ABE:
    ABE_list=ABE_list+Ins_dict[i]

for i in SCNT:
    list1=Ins_dict[i]
    # filter1=list(set(list1) & set(depth40list))
    filter1=list(set(list1)-set(GFP_list))
    filter1=list(set(filter1)-set(ABE_list))
    for j in filter1:
        f2.write('\t'.join([j,'Insertion',i])+'\n')
for i in ABE:
    list1=Ins_dict[i]
    # filter1=list(set(list1) & set(depth40list))
    filter1=list(set(list1)-set(GFP_list))
    filter1=list(set(filter1)-set(SCNT_list))
    for j in filter1:
        f2.write('\t'.join([j,'Insertion',i])+'\n')

####Deletion
GFP_list=[]
for i in GFP:
    GFP_list=GFP_list+del_dict[i]
SCNT_list=[]
for i in SCNT:
    SCNT_list=SCNT_list+del_dict[i]
ABE_list=[]
for i in ABE:
    ABE_list=ABE_list+del_dict[i]

for i in SCNT:
    list1=del_dict[i]
    # filter1=list(set(list1) & set(depth40list))
    filter1=list(set(list1)-set(GFP_list))
    filter1=list(set(filter1)-set(ABE_list))
    for j in filter1:
        f2.write('\t'.join([j,'Deletion',i])+'\n')
for i in ABE:
    list1=del_dict[i]
    # filter1=list(set(list1) & set(depth40list))
    filter1=list(set(list1)-set(GFP_list))
    filter1=list(set(filter1)-set(SCNT_list))
    for j in filter1:
        f2.write('\t'.join([j,'Deletion',i])+'\n')
