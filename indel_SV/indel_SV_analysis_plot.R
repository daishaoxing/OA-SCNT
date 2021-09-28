# venn plot ---------------------------------------------------------------
setwd('/data1/nyy/offtarget/v4indel')
rm(list=ls())
library(venn)
library(reshape2)
library(dplyr)
library(stringr)
library(ggplot2)
library("ggrepel")
##install.packages("ggpubr")
library(ggpubr)
alllist=list()
# i='gatk'
samplelist<-c("090202body","090202head","NC1","NC2","NC3","NC4","NC5","NC6","NC7","NC8","NC9",
              "cloning-JR","cloning-PF","cloning-W",
              "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
for (i in c('gatk','lofreq2','platypus','strelka2')) {
  fname=paste('merge_26',i,'_filter_indel.txt',sep = '')
  myindel<-read.delim(file =fname,header =F, quote = "",sep = "\t",stringsAsFactors=FALSE)
  colnames(myindel)<-c('indel','type','sample')
  indels<-list()
  for (j in samplelist) {
    indel<-myindel$indel[myindel$sample == j]
    indels<-c(indels,list(indel))
  }
  names(indels)<-samplelist
  print(samplelist)

  #######SCNT
  indelSCNT1<-indels[['NC1']]
  indelSCNT2<-indels[['NC2']]
  indelSCNT3<-indels[['NC3']]
  indelSCNT4<-indels[['NC4']]
  indelSCNT5<-indels[['NC5']]
  indelSCNT6<-indels[['NC6']]
  indelSCNT7<-indels[['NC7']]
  indelSCNT8<-indels[['NC8']]
  indelSCNT9<-indels[['NC9']]
  
  indel090202body<-indels[['090202body']]
  indel090202head<-indels[['090202head']]
  
  #######ABE_tissueABE_tissue<-unique(c(snps[['cloning-JR']],snps[['cloning-PF']],snps[['cloning-W']]))
  indelABEJR<-indels[['cloning-JR']]
  indelABEPF<-indels[['cloning-PF']]
  indelABEW<-indels[['cloning-W']]
  
  #######ABE
  indelABE1<-indels[['PC1']]
  indelABE2<-indels[['PC2']]
  indelABE3<-indels[['PC3']]
  indelABE4<-indels[['PC4']]
  indelABE5<-indels[['PC5']]
  indelABE6<-indels[['PC6']]
  indelABE7<-indels[['PC7']]
  indelABE8<-indels[['PC8']]
  
  alllist<-c(alllist,list(indel090202body),list(indel090202head),
             list(indelSCNT1),list(indelSCNT2),list(indelSCNT3),list(indelSCNT4),list(indelSCNT5),
             list(indelSCNT6),list(indelSCNT7),list(indelSCNT8),list(indelSCNT9),
             list(indelABE1),list(indelABE2),list(indelABE3),list(indelABE4),list(indelABE5),
             list(indelABE6),list(indelABE7),list(indelABE8),
             list(indelABEJR),list(indelABEPF),list(indelABEW))
}

gatklist=alllist[1:22]
lofreqlist=alllist[23:44]
platypuslist=alllist[45:66]
strelka2list=alllist[67:88]


SNP_4method=list()
SNPname<-c('090202body','090202head',
           'SCNT1','SCNT2','SCNT3','SCNT4','SCNT5',
           'SCNT6','SCNT7','SCNT8','SCNT9',
           'ABE1','ABE2','ABE3','ABE4','ABE5',
           'ABE6','ABE7','ABE8',
           'ABEJR','ABEPF','ABEW')

for (i in 1:22) {
  gatk=gatklist[[i]]
  Lofreq2=lofreqlist[[i]]
  platypus=platypuslist[[i]]
  strelka2=strelka2list[[i]]
  ######venn
  A1 = gatk
  B1 = Lofreq2
  C1 = platypus
  D1 = strelka2
  intersect1<-Reduce(intersect,list(A1,B1,C1,D1))
  print (length(intersect1))
  A<-paste('gatk(',length(A1),' SNPs)',sep='')
  B<-paste('Lofreq2(',length(B1),' SNPs)',sep='')
  C<-paste('platypus(',length(C1),' SNPs)',sep='')
  D<-paste('strelka2(',length(D1),' SNPs)',sep='')
  sname=paste(A,B,C,D,sep=',')
  x = list(A1,B1,C1,D1)   #将4组数据也就是4个 factor 放入一个 list 中
  outname<-paste('./fig/',SNPname[i],'_DP40_filter_indel_overlap_venn.pdf',sep='')
  pdf(outname, width = 8, height = 8)
  venn(x, snames = sname,ellipse = F, zcolor = "style", cexil = 1, cexsn = 0.8,ilcs=1.5)
  dev.off()
  SNP_4method<-c(SNP_4method,list(intersect1))
}

names(SNP_4method)<-SNPname
df_snp_sample<-data.frame(SNP='SNP',sample='smaple',stringsAsFactors = F)
for (i in 1:length(SNP_4method)) {
  samplesnp<-SNP_4method[[i]]
  samplename<-names(SNP_4method)[i]
  tmpdf<-data.frame(SNP=samplesnp,sample=samplename)
  df_snp_sample<-rbind(df_snp_sample,tmpdf,stringsAsFactors=F)
}

df_snp_sample<-df_snp_sample[-1,]
write.table(df_snp_sample,file = '4method_overlap_DP40_filter_indels.txt',sep = '\t',row.names = F)

indelnumer<-c()
for (i in 1:length(SNP_4method)){
  indelnumer<-c(indelnumer,length(SNP_4method[[i]]))
}
df_indel<-data.frame(number=indelnumer,sample=names(SNP_4method),stringsAsFactors = F)
write.table(df_indel,file = '4method_overlap_DP40_filter_indel_number.txt',sep = '\t',row.names = F)

# ############################################geom_bar --------------------
setwd('/data1/nyy/offtarget/v4indel')
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file = '4method_overlap_DP40_filter_indel_number.txt',sep = '\t',header = T,stringsAsFactors = F)
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
mydata$sample<-factor(mydata$sample, levels=unique(mydata$sample)) 
p<-ggplot(data = mydata, mapping = aes(x = factor(sample), y = number,fill = sample)) + geom_bar(stat = 'identity', position = 'dodge')+
  scale_y_continuous(breaks =breaks)+
  labs(x = "",y = "Number of indels")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_blank())
p
plotfile='./fig/merge_26sample_filter_indel_4method_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)
plotfile='./fig/merge_26sample_filter_indel_4method_count.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)


# #########boxplot --------------------------------------------------------
mydata$Status<-case_when(
  mydata$sample %in% c('ABEJR','ABEPF','ABEW') ~ "SCNT-ABE",
  mydata$sample %in% c("090202body", "090202head",'SCNT1','SCNT2','SCNT3','SCNT4','SCNT5','SCNT6','SCNT7','SCNT8','SCNT9') ~ "SCNT",
  mydata$sample %in% c('ABE1','ABE2','ABE3','ABE4','ABE5','ABE6','ABE7','ABE8') ~ "SCNT-ABE"
)

mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))

############only sum
sem <- function(x) sqrt(var(x)/length(x))
(ABE_mean<-mean(mydata$number[mydata$Status=='SCNT-ABE'])) ##43.72727
(ABE_sem<-sem(mydata$number[mydata$Status=='SCNT-ABE'])) ##9.285464
(SCNT_meam<-mean(mydata$number[mydata$Status=='SCNT'])) ##48.90909
(SCNT_sem<-sem(mydata$number[mydata$Status=='SCNT']))   ##5.649545

(ABE_sd<-sd(mydata$number[mydata$Status=='SCNT-ABE'])) ##30.7964
(SCNT_sd<-sd(mydata$number[mydata$Status=='SCNT']))   ##18.73742
 
var.test(c(mydata$number[mydata$Status=='SCNT-ABE']),c(mydata$number[mydata$Status=='SCNT']))

# F test to compare two variances
# 
# data:  c(mydata$number[mydata$Status == "SCNT-ABE"]) and c(mydata$number[mydata$Status == "SCNT"])
# F = 21.833, num df = 10, denom df = 10, p-value = 0.000035
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   5.87428 81.15049
# sample estimates:
#   ratio of variances 
# 21.83348
shapiro.test(c(mydata$number[mydata$Status=='SCNT-ABE'])) ##p-value = 0.00465
shapiro.test(c(mydata$number[mydata$Status=='SCNT']))     ##p-value = 0.05605


breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
mydata$Status
p<-ggplot(mydata, aes(x=Status, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of indels")+
  scale_y_continuous(breaks = breaks)+
  guides(color = "none", fill = "none")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum-1,hide.ns = TRUE)

p
plotfile='./fig/4method_overlap_filter_indel_sum_boxplot_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_filter_indel_sum_boxplot_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

library("ggrepel")
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center')+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="point")+
  labs(x = "",y = "Number of indels")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  guides(color = "none", fill = "none")+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)

p1<-p+geom_text_repel(data=mydata, aes(label=sample),col="black",alpha = 1)
p
p1

plotfile='./fig/4method_overlap_filter_indel_sum_dotplot_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_filter_indel_sum_dotplot_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

plotfile='./fig/4method_overlap_filter_indel_sum_dotplot_wilcoxlabel.pdf'
ggsave(plotfile, plot=p1, dpi = 600,width = 10, height = 6)
plotfile='./fig/4method_overlap_filter_indel_sum_dotplot_wilcoxlabel.png'
ggsave(plotfile, plot=p1, dpi = 600,width = 10, height = 6)


# ######delly SV  ---------------------------------------------------------
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library("ggrepel")
svdata1<-read.delim(file = '/data1/nyy/offtarget_addfq/delly_SV/merged_SV_filter_count.tab',sep = '\t',header = F,stringsAsFactors = F)
colnames(svdata1)<-c('pos','type','sample')
samplelist<-unique(svdata1$sample)

svdata<-svdata1 %>%
  dplyr::group_by(type,sample) %>%
  dplyr::summarise(number = n())

svdata$Status<-case_when(
  svdata$sample %in% c("cloning-JR","cloning-PF","cloning-W") ~ "SCNT-ABE",
  svdata$sample %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  svdata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
svdata$Status<-factor(svdata$Status, levels=c('SCNT-ABE', 'SCNT'))
svdata$type<-factor(svdata$type, levels=c("DEL","INS","DUP","INV","BND"))

####all sv
mydata<-svdata
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=type, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  scale_x_discrete(labels=c("Deletion","Insertion","Duplication","Inversion","Translocation"))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
# stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
#                    method = "wilcox.test",method.args =list(alternative = "greater"),size=5,label.y=maximum-1,hide.ns = TRUE)
p
plotfile='./fig/delly_SV_boxplot.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)
plotfile='./fig/delly_SV_boxplot.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)

##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
####Deletion
mydata<-svdata[svdata$type=='DEL',]
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center')+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="point")+
  labs(x = "",y = "Number of Deletion")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=16))+
  guides(color = "none", fill = "none")+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/delly_SV_del_dotplot.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)
plotfile='./fig/delly_SV_del_dotplot.png'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)

####Insertion
mydata<-svdata[svdata$type=='INS',]
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center')+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="point")+
  labs(x = "",y = "Number of Insertion")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=16))+
  guides(color = "none", fill = "none")+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/delly_SV_INS_dotplot.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)
plotfile='./fig/delly_SV_INS_dotplot.png'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)

####Duplication
mydata<-svdata[svdata$type=='DUP',]
mydata$number<-mydata$number/1000
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center')+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="point")+
  labs(x = "",y = "Number of Duplication (X1000)")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=16))+
  guides(color = "none", fill = "none")+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/delly_SV_DUP_dotplot.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)
plotfile='./fig/delly_SV_DUP_dotplot.png'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)

####Inversion
mydata<-svdata[svdata$type=='INV',]
mydata$number<-mydata$number/1000
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center')+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="point")+
  labs(x = "",y = "Number of Inversion (X1000)")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=16))+
  guides(color = "none", fill = "none")+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/delly_SV_INV_dotplot.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)
plotfile='./fig/delly_SV_INV_dotplot.png'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)

####Translocation
mydata<-svdata[svdata$type=='BND',]
breaks<-pretty(range(0,mydata$number), 6)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center')+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="point")+
  labs(x = "",y = "Number of Translocation")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=16))+
  guides(color = "none", fill = "none")+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/delly_SV_BND_dotplot.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)
plotfile='./fig/delly_SV_BND_dotplot.png'
ggsave(plotfile, plot=p, dpi = 600,width = 4, height = 4)
