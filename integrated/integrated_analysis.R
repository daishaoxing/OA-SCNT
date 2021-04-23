
##############
rm(list=ls())
library(ggplot2)
# install.packages('venn')
library(venn)
library(reshape2)
setwd("/home/devdata/nyy/nyy_GFP/GFP/integrated")
DNAdata<-read.table(file ='DNA_A2Ggene.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
RNAdata<-read.table(file ='RNA_A2Ggene.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(RNAdata)<-c('gene','sample')
colnames(DNAdata)<-c('gene','sample')
NClist<-c('NC1','NC2','NC3','NC4','NC5')
DNAdata<-DNAdata[! DNAdata$sample %in% NClist,]
RNAdata<-RNAdata[! RNAdata$sample %in% NClist,]
unique(DNAdata$sample)
unique(RNAdata$sample)
# [1] "PC1"       "PC2"       "PC3"       "PC4"       "PC5"       "PC6"       "PC7"       "PC8"      
# [9] "cloningW"  "cloningPF" "cloningJR"

DNA=unique(DNAdata$gene)
RNA=unique(RNAdata$gene)
intersect1<-Reduce(intersect,list(DNA,RNA))
print (length(intersect1))

# print (c(length(A1),length(B1),length(C1),length(D1)))
A<-paste('DNA_edit (',length(DNA),' genes)',sep='')
B<-paste('RNA_edit (',length(RNA),' genes)',sep='')

sname=paste(A,B,sep=',')
x = list(DNA,RNA)   #将2组数据也就是2个 factor 放入一个 list 中
pdf(file = "DNA_RNA_gene_overlap_venn.pdf",
    width = 4,
    height = 4)
venn(x, snames = sname, zcolor = "style", cexil = 1,ilcs=1)
dev.off()

##############site overlap##############
rm(list=ls())
library(ggplot2)
# install.packages('venn')
library(venn)
library(reshape2)
setwd("/home/devdata/nyy/nyy_GFP/GFP/integrated")
DNAdata<-read.table(file ='DNA_A2Gsite.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
RNAdata<-read.table(file ='RNA_A2Gsite.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(RNAdata)<-c('site','sample')
colnames(DNAdata)<-c('site','sample')
NClist<-c('NC1','NC2','NC3','NC4','NC5')
DNAdata<-DNAdata[! DNAdata$sample %in% NClist,]
RNAdata<-RNAdata[! RNAdata$sample %in% NClist,]
unique(DNAdata$sample)
unique(RNAdata$sample)

DNA=unique(DNAdata$site)
RNA=unique(RNAdata$site)
intersect1<-Reduce(intersect,list(DNA,RNA))
print (length(intersect1))

# print (c(length(A1),length(B1),length(C1),length(D1)))
A<-paste('DNA sites (',length(DNA),' sites)',sep='')
B<-paste('RNA sites (',length(RNA),' sites)',sep='')

sname=paste(A,B,sep=',')
x = list(DNA,RNA)   #将5组数据也就是5个 factor 放入一个 list 中
pdf(file = "DNA_RNA_site_overlap_venn.pdf",
    width = 6,
    height = 4)
venn(x, snames = sname, zcolor = "style", cexil = 1,ilcs=1)
dev.off()

