import os
#############background 81 macaca DP>20
f1=open('/nas01/nyy/vcffrom27/population_all_SNP.txt','r')
population_pos=[]
for i in f1.readlines():
    tmp=i.strip().split('\t')
    population_pos.append(tmp[0])
f1.close()

f1=open('/nas01/nyy/vcffrom27/population_all_SNP_p2.txt','r')
for i in f1.readlines():
    tmp=i.strip().split('\t')
    tmp=i.strip().split('\t')
    population_pos.append(tmp[0])
f1.close()

f1=open('/nas01/nyy/vcffrom27/population_all_SNP_p3.txt','r')
for i in f1.readlines():
    tmp=i.strip().split('\t')
    tmp=i.strip().split('\t')
    population_pos.append(tmp[0])
f1.close()
print(population_pos[3000])


path='/data1/nyy/offtarget/vcf_2021'
method=['gatk','lofreq2','platypus','strelka2']
sample=['08431GFP','08431neg','08431pos','08431WT','090202body','090202head','cloning-JR','cloning-PF','cloning-W','NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8']
dp=20
f2=open('GATK_DP20_SNP_count.txt','w')
f2.write('\t'.join(['Individual','nofilter','filter'])+'\n')
dict1={}
for i in sample:
    fname=os.path.join(path,i+'_gatk_SNP.txt')
    print(fname)
    f1=open(fname,'r')
    for j in f1.readlines()[1:]:###POS     REF     ALT     GT      AF      AF2     DP      info 
        tmp=j.strip().split('\t')
        if float(tmp[6])>=dp:
            str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
            dict1.setdefault(i,[]).append(str1)
    f1.close()
dict2={}
for i in dict1.keys():
    filterSNP=set(dict1[i])-set(population_pos)
    number1=len(dict1[i])
    number2=len(filterSNP)
    dict2[i]=[number1,number2]
    f2.write('\t'.join([i,str(number1),str(number2)])+'\n')
print(dict2)
f2.close()

dp=40
f2=open('GATK_DP40_SNP_count.txt','w')
f2.write('\t'.join(['Individual','nofilter','filter'])+'\n')
dict1={}
for i in sample:
    fname=os.path.join(path,i+'_gatk_SNP.txt')
    print(fname)
    f1=open(fname,'r')
    for j in f1.readlines()[1:]:###POS     REF     ALT     GT      AF      AF2     DP      info 
        tmp=j.strip().split('\t')
        if float(tmp[6])>=dp:
            str1=tmp[0]+'|'+tmp[1]+'>'+tmp[2]
            dict1.setdefault(i,[]).append(str1)
    f1.close()
dict2={}
for i in dict1.keys():
    filterSNP=set(dict1[i])-set(population_pos)
    number1=len(dict1[i])
    number2=len(filterSNP)
    dict2[i]=[number1,number2]
    f2.write('\t'.join([i,str(number1),str(number2)])+'\n')
print(dict2)
f2.close()