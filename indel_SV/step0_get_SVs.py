f1=open('merged_DEL_new.vcf','r')
f2=open('merged_DEL_new_info.tab','w')
f2.write("chrpos\t08431GFP\t08431neg\t08431pos\t08431WT\t090202body\tNC1\tNC2\tNC3\tNC4\tNC5\tNC6\tNC7\tNC8\tNC9\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\t090202head\tcloning-JR\tcloning-PF\tcloning-W\n")
for i in f1.readlines():
    if i[0] !='#':
        tmp=i.strip().split('\t')
        chrpos=tmp[0]+':'+tmp[1]
        QUAL=tmp[5]
        FILTER=tmp[6]
        INFO=tmp[7].split(';')
        PRECISE=INFO[0]
        if float(QUAL)>20 and FILTER=='PASS' and PRECISE=='PRECISE':
            info_dict={}
            for k in INFO[1:-1]:
                (a,b)=k.split('=')
                info_dict[a]=b
            PE=info_dict['PE']
            SR=info_dict['SR']
            MAPQ=info_dict['MAPQ']
            if float(SR) >0 or float(PE)>3:
                list1=[]
                for j in tmp[9:]:
                    FORMAT=j.strip().split(':')
                    gt='0'
                    GQ=FORMAT[2]
                    FT=FORMAT[3]
                    if '1' in FORMAT[0]:
                        gt='1'
                    list1.append(gt)
                f2.write(chrpos+'\t'+'\t'.join(list1)+'\n')
        else:
            if INFO[-1] =='IMPRECISE':
                print(INFO)
f1.close()
f2.close()


f1=open('merged_DUP_new.vcf','r')
f2=open('merged_DUP_new_info.tab','w')
f2.write("chrpos\t08431GFP\t08431neg\t08431pos\t08431WT\t090202body\tNC1\tNC2\tNC3\tNC4\tNC5\tNC6\tNC7\tNC8\tNC9\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\t090202head\tcloning-JR\tcloning-PF\tcloning-W\n")
for i in f1.readlines():
    if i[0] !='#':
        tmp=i.strip().split('\t')
        chrpos=tmp[0]+':'+tmp[1]
        QUAL=tmp[5]
        FILTER=tmp[6]
        INFO=tmp[7].split(';')
        PRECISE=INFO[0]
        if float(QUAL)>20 and FILTER=='PASS' and PRECISE=='PRECISE':
            info_dict={}
            for k in INFO[1:-1]:
                (a,b)=k.split('=')
                info_dict[a]=b
            PE=info_dict['PE']
            SR=info_dict['SR']
            MAPQ=info_dict['MAPQ']
            if float(SR) >0 or float(PE)>3:
                list1=[]
                for j in tmp[9:]:
                    FORMAT=j.strip().split(':')
                    gt='0'
                    GQ=FORMAT[2]
                    FT=FORMAT[3]
                    if '1' in FORMAT[0]:
                        gt='1'
                    list1.append(gt)
                f2.write(chrpos+'\t'+'\t'.join(list1)+'\n')
        else:
            if INFO[-1] =='IMPRECISE':
                print(INFO)
f1.close()
f2.close()


f1=open('merged_INV_new.vcf','r')
f2=open('merged_INV_new_info.tab','w')
f2.write("chrpos\t08431GFP\t08431neg\t08431pos\t08431WT\t090202body\tNC1\tNC2\tNC3\tNC4\tNC5\tNC6\tNC7\tNC8\tNC9\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\t090202head\tcloning-JR\tcloning-PF\tcloning-W\n")
for i in f1.readlines():
    if i[0] !='#':
        tmp=i.strip().split('\t')
        chrpos=tmp[0]+':'+tmp[1]
        QUAL=tmp[5]
        FILTER=tmp[6]
        INFO=tmp[7].split(';')
        PRECISE=INFO[0]
        if float(QUAL)>20 and FILTER=='PASS' and PRECISE=='PRECISE':
            info_dict={}
            for k in INFO[1:-1]:
                (a,b)=k.split('=')
                info_dict[a]=b
            PE=info_dict['PE']
            SR=info_dict['SR']
            MAPQ=info_dict['MAPQ']
            if float(SR) >0 or float(PE)>3:
                list1=[]
                for j in tmp[9:]:
                    FORMAT=j.strip().split(':')
                    gt='0'
                    GQ=FORMAT[2]
                    FT=FORMAT[3]
                    if '1' in FORMAT[0]:
                        gt='1'
                    list1.append(gt)
                f2.write(chrpos+'\t'+'\t'.join(list1)+'\n')
        else:
            if INFO[-1] =='IMPRECISE':
                print(INFO)
f1.close()
f2.close()

f1=open('merged_INS_new.vcf','r')
f2=open('merged_INS_new_info.tab','w')
f2.write("chrpos\t08431GFP\t08431neg\t08431pos\t08431WT\t090202body\tNC1\tNC2\tNC3\tNC4\tNC5\tNC6\tNC7\tNC8\tNC9\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\t090202head\tcloning-JR\tcloning-PF\tcloning-W\n")
for i in f1.readlines():
    if i[0] !='#':
        tmp=i.strip().split('\t')
        chrpos=tmp[0]+':'+tmp[1]
        QUAL=tmp[5]
        FILTER=tmp[6]
        INFO=tmp[7].split(';')
        PRECISE=INFO[0]
        if float(QUAL)>20 and FILTER=='PASS' and PRECISE=='PRECISE':
            info_dict={}
            for k in INFO[1:-1]:
                (a,b)=k.split('=')
                info_dict[a]=b
            PE=info_dict['PE']
            SR=info_dict['SR']
            MAPQ=info_dict['MAPQ']
            if float(SR) >0 or float(PE)>3:
                list1=[]
                for j in tmp[9:]:
                    FORMAT=j.strip().split(':')
                    gt='0'
                    GQ=FORMAT[2]
                    FT=FORMAT[3]
                    if '1' in FORMAT[0]:
                        gt='1'
                    list1.append(gt)
                f2.write(chrpos+'\t'+'\t'.join(list1)+'\n')
        else:
            if INFO[-1] =='IMPRECISE':
                print(INFO)
f1.close()
f2.close()


f1=open('merged_BND_new.vcf','r')
f2=open('merged_BND_new_info.tab','w')
f2.write("chrpos\t08431GFP\t08431neg\t08431pos\t08431WT\t090202body\tNC1\tNC2\tNC3\tNC4\tNC5\tNC6\tNC7\tNC8\tNC9\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\t090202head\tcloning-JR\tcloning-PF\tcloning-W\n")

for i in f1.readlines():
    if i[0] !='#':
        tmp=i.strip().split('\t')
        chrpos=tmp[0]+':'+tmp[1]
        QUAL=tmp[5]
        FILTER=tmp[6]
        INFO=tmp[7].split(';')
        PRECISE=INFO[0]
        if float(QUAL)>20 and FILTER=='PASS' and PRECISE=='PRECISE':
            info_dict={}
            for k in INFO[1:-1]:
                (a,b)=k.split('=')
                info_dict[a]=b
            PE=info_dict['PE']
            SR=info_dict['SR']
            MAPQ=info_dict['MAPQ']
            if float(SR) >0 or float(PE)>3:
                list1=[]
                for j in tmp[9:]:
                    FORMAT=j.strip().split(':')
                    gt='0'
                    GQ=FORMAT[2]
                    FT=FORMAT[3]
                    if '1' in FORMAT[0]:
                        gt='1'
                    list1.append(gt)
                f2.write(chrpos+'\t'+'\t'.join(list1)+'\n')
        else:
            if INFO[-1] =='IMPRECISE':
                print(INFO)
f1.close()
f2.close()


####CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	08431GFP	08431neg	08431pos	08431WT	090202body	NC1	NC2	NC3	NC4	NC5	NC6	NC7	NC8	NC9	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	090202head	cloning-JR	cloning-PF	cloning-W
###1       128527542       DEL00000981     T       <DEL>   .       PASS    PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.7.6;CHR2=1;END=128527847;PE=2;MAPQ=60;CT=3to5;CIPOS=-23,23;CIEND=-23,23;INSLEN=0;HOMLEN=22;SR=6;SRQ=0.994413;CONSENSUS=TTACAGGAATGAGCCACTGTGCCCGGCCCCAGCCAACCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAAGACGGAGTCTCGCTCTGTTGCCCAAGTTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAACCTCTGCCTCTGGGGTTCAAGTGATTCTCCTGCCTCAGCCTCCTAATTA;CE=1.93772      GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV     ./.:.:.:.:.:.:.:.:.:.:.:.       1/1:-13.299,-1.20311,0:12:LowQual:55:0:42:0:6:3:0:4     ./.:.:.:.:.:.:.:.:.:.:.:../.:.:.:.:.:.:.:.:.:.:.:.        ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:. ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:. ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:. ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:../.:.:.:.:.:.:.:.:.:.:.:.        ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.       ./.:.:.:.:.:.:.:.:.:.:.:.
# SVs are flagged as FILTER:LowQual if PE <3 OR MAPQ <20 (for translocations: PE <5 OR MAPQ <20), otherwise, the SV results in a FILTER:PASS. PRECISE variants will have split-read support (SR >0).
# QUAL >=20; SR >0; DP>3; marked PASS and PRECISE tag; Bing Su study
# Site Filtering QUAL >=20; (SR >0 or PE>3) MAPQ >20 marked PASS and PRECISE tag; 
# Genotype Filtering GT:GL:GQ:FT:   GQ>30  FT==PASS