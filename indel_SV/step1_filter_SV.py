import os
import _pickle as cPickle
flist=['DEL','DUP','INV','INS','BND'] ##merged_INS_new_info.tab
# flist=['DEL','DUP','INV','BND'] ##merged_INS_new_info.tab

f2=open('merged_SV_filter_count.tab','w')
for k in flist:
    fname='merged_'+k+'_new_info.tab'
    f1=open(fname,'r')
    print(k)
    # depth40list = cPickle.load(open("/nas01/nyy/vcffrom27/depth40_vcfpos.pkl","rb"))
    alllist=['08431GFP', '08431neg', '08431pos', '08431WT', '090202body', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5', 'NC6', 'NC7', 'NC8', 'NC9', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', '090202head', 'cloning-JR', 'cloning-PF', 'cloning-W']
    GFP=['08431neg','08431pos','08431GFP','08431WT']
    SCNT=['090202body','090202head','NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9']
    ABE=['cloning-JR','cloning-PF','cloning-W','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8']
    del_dict={}
    for i in f1.readlines()[1:]:
        tmp=i.strip().split('\t')
        list1=tmp[1:]
        for j in range(len(list1)):
            if list1[j]=='1':
                del_dict.setdefault(alllist[j],[]).append(tmp[0])
    # print(del_dict['08431neg'])
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
            f2.write('\t'.join([j,k,i])+'\n')
    for i in ABE:
        list1=del_dict[i]
        # filter1=list(set(list1) & set(depth40list))
        filter1=list(set(list1)-set(GFP_list))
        filter1=list(set(filter1)-set(SCNT_list))
        for j in filter1:
            f2.write('\t'.join([j,k,i])+'\n')
    f1.close()
f2.close()
