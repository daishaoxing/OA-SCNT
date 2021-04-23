# install.packages('venn')
setwd('/home/devdata/nyy/nyy_GFP/GFP/filter_vcf')


# Reduce(intersect,  list(v1 = c("a","b","c","d"),v2 = c("a","b","e"),
#                         v3 = c("a","f","g"),v4 = c("a","f","e","g"))
#      )


########
rm(list=ls())
library(ggplot2)
library(venn)
library(reshape2)
#options(scipen = 9)
mySNP<-read.table(file ='merge_15sample_filter_add_4method_overlap_venn.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mySNP)<-c('sample','SNP','method')
# unique(mySNP$sample)
# [1] "cloningPF" "PC5"       "NC1"       "PC4"       "PC1"       "PC7"       "PC6"       "cloningJR"
# [9] "NC2"       "NC3"       "PC3"       "PC8"       "cloningW"  "PC2" 
samplelist<-c('NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8',"cloningPF","cloningJR","cloningW")
# samplelist<-c('NC1')
for (i in samplelist){
  A1 = mySNP$SNP[mySNP$method == 'lofreq2' & mySNP$sample == i] #把同 sample 的 gene 们归类到一个 list中，下略
  B1 = mySNP$SNP[mySNP$method == 'platypus' & mySNP$sample == i]
  C1 = mySNP$SNP[mySNP$method == 'gatk' & mySNP$sample == i]
  D1 = mySNP$SNP[mySNP$method == 'strelka' & mySNP$sample == i]
  intersect1<-Reduce(intersect,list(A1,B1,C1,D1))
  # print (intersect1)
  print (length(intersect1))
  # print (c(length(A1),length(B1),length(C1),length(D1)))
  A<-paste('Lofreq2(',length(A1),' SNPs)',sep='')
  B<-paste('Platypus(',length(B1),' SNPs)',sep='')
  C<-paste('gatk(',length(C1),' SNPs)',sep='')
  D<-paste('Strelka2(',length(D1),' SNPs)',sep='')
  sname=paste(A,B,C,D,sep=',')
  x = list(A1,B1,C1,D1)   #将4组数据也就是4个 factor 放入一个 list 中
  outname<-paste('./fig/',i,'_SNP_overlap_venn.tif',sep='')
  tiff(outname, res = 300,width = 10, height = 8, units = 'in')
  venn(x, snames = sname,ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8,ilcs=2)
  dev.off()
}


#######
rm(list=ls())
library(ggplot2)
library(venn)
library(reshape2)
#options(scipen = 9)
samplelist<-c('07027', '08431GFP', '08431WT', 'NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','cloning-W','cloning-PF','cloning-JR')
for (i in samplelist){
  snpfile<-paste(i,'_SNP_overlap_venn.txt',sep='')
  mySNP<-read.table(file =snpfile,header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  colnames(mySNP)<-c('SNP','method')
  A1 = mySNP$SNP[mySNP$method == 'lofreq2'] #把同 sample 的 gene 们归类到一个 list中，下略
  B1 = mySNP$SNP[mySNP$method == 'platypus']
  C1 = mySNP$SNP[mySNP$method == 'gatk']
  D1 = mySNP$SNP[mySNP$method == 'strelka2']
  c(length(A1),length(B1),length(C1),length(D1))
  A<-paste('Lofreq2(',length(A1),' SNPs)',sep='')
  B<-paste('Platypus(',length(B1),' SNPs)',sep='')
  C<-paste('gatk(',length(C1),' SNPs)',sep='')
  D<-paste('Strelka2(',length(D1),' SNPs)',sep='')
  sname=paste(A,B,C,D,sep=',')
  x = list(A1,B1,C1,D1)   #将4组数据也就是4个 factor 放入一个 list 中
  outname<-paste('./fig/',i,'_SNP_overlap_venn.tif',sep='')
  tiff(outname, res = 400,width = 10, height = 10, units = 'in')
  venn(x, snames = sname, ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8)
  dev.off()
}

####find ./ -regex '.*\(jpg\|JPG\|png\|PNG\|tif\|jpeg\)' -size +50k -exec convert -resize 50%x50% {} {} \;
####find ./ -iname '.*filter_overlap_venn.tif' -exec convert -resize "50%x50%" {} {} \;

rm(list=ls())
library(ggplot2)
library(venn)
library(reshape2)
#options(scipen = 9)
samplelist<-c('NC1', 'NC2', 'NC3', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8')
for (i in samplelist){
  snpfile<-paste(i,'_SNP_filter_overlap_venn.txt',sep='')
  mySNP<-read.table(file =snpfile,header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  colnames(mySNP)<-c('SNP1','SPN2','method')
  mySNP$SNP<-paste(mySNP$SNP1,mySNP$SNP2,sep='|')
  A1 = mySNP$SNP[mySNP$method == 'lofreq2'] #把同 sample 的 gene 们归类到一个 list中，下略
  B1 = mySNP$SNP[mySNP$method == 'platypus']
  C1 = mySNP$SNP[mySNP$method == 'freebayes']
  D1 = mySNP$SNP[mySNP$method == 'strelka2']
  c(length(A1),length(B1),length(C1),length(D1))
  A<-paste('Lofreq2(',length(A1),' SNPs)',sep='')
  B<-paste('Platypus(',length(B1),' SNPs)',sep='')
  C<-paste('Freebayes(',length(C1),' SNPs)',sep='')
  D<-paste('Strelka2(',length(D1),' SNPs)',sep='')
  sname=paste(A,B,C,D,sep=',')
  x = list(A1,B1,C1,D1)   #将4组数据也就是4个 factor 放入一个 list 中
  outname<-paste('./fig/',i,'_SNP_filter_overlap_venn.tif',sep='')
  tiff(outname, res = 400,width = 10, height = 10, units = 'in')
  venn(x, snames = sname, ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8)
  dev.off()
}

##########
rm(list=ls())
library(ggplot2)
library(venn)
library(reshape2)
#options(scipen = 9)
mySNP<-read.table(file ='sample17_overlap_SNP_venn.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mySNP)<-c('SNP','sample')
NC1 = mySNP$SNP[mySNP$sample == 'NC1']
NC2 = mySNP$SNP[mySNP$sample == 'NC2']
NC3 = mySNP$SNP[mySNP$sample == 'NC3']
WT = mySNP$SNP[mySNP$sample == '08431GFP']
c(length(NC1),length(NC2),length(NC3),length(WT))
A<-paste('NC1(',length(NC1),' SNPs)',sep='')
B<-paste('NC2(',length(NC2),' SNPs)',sep='')
C<-paste('NC3(',length(NC3),' SNPs)',sep='')
D<-paste('GFP08431(',length(WT),' SNPs)',sep='')
sname=paste(A,B,C,D,sep=',')
x = list(NC1,NC2,NC3,WT)   #将4组数据也就是4个 factor 放入一个 list 中
outname<-paste('./fig/','GFP_vs_NC_SNP_venn.tif',sep='')
tiff(outname, res = 400,width = 10, height = 10, units = 'in')
venn(x, snames = sname, ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8)
dev.off()


PC1 = mySNP$SNP[mySNP$sample == 'PC1']
PC2 = mySNP$SNP[mySNP$sample == 'PC2']
PC3 = mySNP$SNP[mySNP$sample == 'PC3']
PC4 = mySNP$SNP[mySNP$sample == 'PC4']
WT = mySNP$SNP[mySNP$sample == '08431GFP']
c(length(PC1),length(PC2),length(PC3),length(PC4),length(WT))
A<-paste('PC1(',length(PC1),' SNPs)',sep='')
B<-paste('PC2(',length(PC2),' SNPs)',sep='')
C<-paste('PC3(',length(PC3),' SNPs)',sep='')
D<-paste('PC4(',length(PC4),' SNPs)',sep='')
E<-paste('GFP08431(',length(WT),' SNPs)',sep='')
sname=paste(A,B,C,D,E,sep=',')
x = list(PC1,PC2,PC3,PC4,WT)   #将4组数据也就是4个 factor 放入一个 list 中
outname<-paste('./fig/','GFP_vs_PC14_SNP_venn.tif',sep='')
tiff(outname, res = 400,width = 10, height = 10, units = 'in')
venn(x, snames = sname, ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8)
dev.off()

##############
PC5 = mySNP$SNP[mySNP$sample == 'PC5']
PC6 = mySNP$SNP[mySNP$sample == 'PC6']
PC7 = mySNP$SNP[mySNP$sample == 'PC7']
PC8 = mySNP$SNP[mySNP$sample == 'PC8']
WT = mySNP$SNP[mySNP$sample == '08431GFP']
c(length(PC5),length(PC6),length(PC7),length(PC8),length(WT))
A<-paste('PC5(',length(PC5),' SNPs)',sep='')
B<-paste('PC6(',length(PC6),' SNPs)',sep='')
C<-paste('PC7(',length(PC7),' SNPs)',sep='')
D<-paste('PC8(',length(PC8),' SNPs)',sep='')
E<-paste('GFP08431(',length(WT),' SNPs)',sep='')
sname=paste(A,B,C,D,E,sep=',')
x = list(PC5,PC6,PC7,PC8,WT)   #将4组数据也就是4个 factor 放入一个 list 中
outname<-paste('./fig/','GFP_vs_PC58_SNP_venn.tif',sep='')
tiff(outname, res = 400,width = 10, height = 10, units = 'in')
venn(x, snames = sname, ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8)
dev.off()

#############all NC, PC
allNC<-c(NC1,NC2,NC3)
allPC<-c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8)
allNC<-unique(allNC)
allPC<-unique(allPC)
WT = mySNP$SNP[mySNP$sample == '08431GFP']
c(length(allNC),length(allPC),length(WT))
A<-paste('allNC(',length(allNC),'SNPs)',sep='')
B<-paste('allPC(',length(allPC),'SNPs)',sep='')
C<-paste('GFP08431(',length(WT),'SNPs)',sep='')
sname=paste(A,B,C,sep=',')
x = list(allNC,allPC,WT)   #将3组数据也就是3个 factor 放入一个 list 中
outname<-paste('./fig/','GFP_vs_all_NC_PC_SNP_venn.tif',sep='')
tiff(outname, res = 400,width = 10, height = 10, units = 'in')
venn(x, snames = sname, zcolor = "style", cexil = 1, cexsn = 0.8)
dev.off()

############PC1-8,NC1
samplelist<-c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8')
for (i in samplelist){
  s <- mySNP$SNP[mySNP$sample == i]
  print (i)
  A<-paste(i,'(',length(s),'SNPs)',sep='')
  B<-paste('NC1(',length(NC1),'SNPs)',sep='')
  C<-paste('GFP08431(',length(WT),'SNPs)',sep='')
  sname=paste(A,B,C,sep=',')
  print (sname)
  x = list(s,NC1,WT)   #将3组数据也就是3个 factor 放入一个 list 中
  outname<-paste('./fig/',i,'_GFP_vs_NC_venn.tif',sep='')
  tiff(outname, res = 400,width = 10, height = 10, units = 'in')
  venn(x, snames = sname, zcolor = "style", cexil = 1, cexsn = 0.8)
  dev.off()
}

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample17_venn_count.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(0,maximum+maximum/25), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('WT07027','WT08431','GFP08431','cloning','NC','PC'))
my_comparisons <- list(c('cloning', 'NC'),c('NC', 'PC'),c('cloning', 'PC'))

p<-ggbarplot(mydata, x = "Status", y = "number",size=1.5,color = "Status",
             add = c("jitter",'mean_se'),merge=TRUE)+
  labs(x = "",y = "Number of SNPs")+
  scale_y_continuous(limits = c(-1,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)
plotfile='./fig/barplot_sample17_SNP_count.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/barplot_sample17_SNP_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlap.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(0,maximum+maximum/40), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggbarplot(mydata, x = "Status", y = "number",size=1.5,color = "Status",
             add = c("jitter",'mean_se'),merge=TRUE)+
  labs(x = "",y = "Number of SNPs")+
  scale_y_continuous(limits = c(-1,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)
plotfile='./fig/barplot_sample11_PC_Vs_NC_SNP_filter.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/barplot_sample11_PC_Vs_NC_SNP_filter.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlapV1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(0,maximum+maximum/40), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggbarplot(mydata, x = "Status", y = "number",size=1.5,color = "Status",
             add = c("jitter",'mean_se'),merge=TRUE)+
  labs(x = "",y = "Number of SNPs")+
  scale_y_continuous(limits = c(-1,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)
plotfile='./fig/barplot_sample11_PC_Vs_NC_SNP_filterV1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/barplot_sample11_PC_Vs_NC_SNP_filterV1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlapV2.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(0,maximum+maximum/40), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggbarplot(mydata, x = "Status", y = "number",size=1.5,color = "Status",
             add = c("jitter",'mean_se'),merge=TRUE)+
  labs(x = "",y = "Number of SNPs")+
  scale_y_continuous(limits = c(-1,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)
plotfile='./fig/barplot_sample11_PC_Vs_NC_SNP_filterV2.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/barplot_sample11_PC_Vs_NC_SNP_filterV2.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

##########barplot geom_boxplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample17_venn_count.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(minimum,maximum+maximum/30), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('WT07027','WT08431','GFP08431','cloning','NC','PC'))
my_comparisons <- list(c('cloning', 'NC'),c('NC', 'PC'),c('cloning', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample17_SNP_count.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample17_SNP_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='all_sample_4method_SNP_overlap_filter20.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[mydata$sampletype %in% c('NC','PC'),c(6,1,7)]
colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(minimum,maximum+maximum/30), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$sampletype, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample17_SNP_count.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample17_SNP_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlap.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(minimum,maximum+maximum/40), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))
p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filter.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filter.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlapV1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
# breaks<-pretty(range(minimum,maximum+maximum/80), 8)
# maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlapV2.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
# breaks<-pretty(range(0,maximum+maximum/40), 8)
# maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV2.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV2.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlapa.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(minimum,maximum+maximum/40), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))
p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filtera.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filtera.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlapV1a.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
# breaks<-pretty(range(minimum,maximum+maximum/80), 8)
# maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV1a.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV1a.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlapV2a.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
breaks<-pretty(range(minimum,maximum+maximum/40), 8)
maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV2a.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_filterV2a.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

##############GFP second filter
####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlap_gfpfilter.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
# breaks<-pretty(range(0,maximum+maximum/40), 8)
# maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_gfpfilter.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_gfpfilter.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

##############GFP second filter_v1
####barplot
rm(list=ls())
library(reshape2)
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata<-read.delim(file ='sample11_SNPs_count_overlap_gfpfilter1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('number','Individual','Status')
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
minimum<- breaks[1]
# breaks<-pretty(range(0,maximum+maximum/40), 8)
# maximum<- breaks[length(breaks)]
mydata$Status<-factor(mydata$Status, levels=c('NC','PC'))
my_comparisons <- list(c('NC', 'PC'))

p<-ggplot(mydata,aes(x=Status, y=number, fill=Status)) +
  geom_boxplot(outlier.colour = NA)+
  #geom_jitter(color="black", size=0.4,position=position_dodge(1)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "Individual",y = "Number of SNPs")+
  scale_y_continuous(limits = c(minimum,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=45)
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=5)

plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_gfpfilter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/geom_boxplot_sample11_PC_Vs_NC_SNP_gfpfilter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

