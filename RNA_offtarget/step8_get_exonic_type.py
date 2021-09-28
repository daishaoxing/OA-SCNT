file_list=['PC1','PC2','PC3','PC4','PC5']
fin=open('merge_GATK_10sample_Jitter.txt','r')
dict1={}
for i in fin.readlines()[1:]:
    tmp=i.strip().split('\t')
    dict1.setdefault(tmp[0],[]).append([tmp[-1],tmp[1]])

for i in file_list:
    fname=i+'_RNA_filter_snps.exonic_variant_function'
    f1=open(fname,'r')
    f2=open(i+'_exonic_variant_infowithtype.txt','w')
    dict2={}
    for j in dict1[i]:
        dict2[j[0]]=j[1]
    # print(dict2)
    for j in f1.readlines():
        tmp=j.strip().split('\t')
        info=tmp[2].split(',')[0]
        (gene,trans,exon,ntalt,aaalt)=info.split(':')
        gene=gene[:-2]
        comments=tmp[-1].split('|') ##comments:|20.369999999999997|11:848308|A>G
        pos_alt=comments[2]+'|'+comments[3]
        if pos_alt in dict2.keys():
            f2.write('\t'.join([tmp[1],gene,exon,pos_alt,aaalt,dict2[pos_alt]])+'\n')
