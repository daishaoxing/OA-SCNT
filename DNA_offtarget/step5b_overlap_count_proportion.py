f1=open('merge_15sample_filter_add_4method_count_merge.txt','r')
f2=open('merge_15sample_filter_add_4method_count_merge_proportion.txt','w')
f2.write("Individual\tbase_change\tnumber\n")

def sumlist(a):
    b=0
    for i in a:
        tmp=i.strip().split('\t')
        b=b+int(tmp[-1])
    return int(b)

lines=f1.readlines()[1:]
sum={}
for i in range(0,84,6):
    b=sumlist(lines[i:(i+6)])
    sum[i]=b
    print (i)
for i in range(0,84,6):
    for j in lines[i:(i+6)]:
        tmp=j.strip().split('\t')
        pro=float(tmp[-1])*100/sum[i]
        print(i,j,sum[i],pro)
        f2.write(tmp[0]+'\t'+tmp[1]+'\t'+str(round(pro,2))+'\n')
f1.close()
f2.close()

