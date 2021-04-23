#######new PCR version
rm(list=ls())
# installed_info<-installed.packages()
# install.packages("ggpubr",lib="/home/dsx/R/x86_64-pc-linux-gnu-library/4.0")
# remove.packages('ggpubr',lib = '/home/dsx/R/x86_64-pc-linux-gnu-library/4.0')
library(ggplot2)
library(ggpubr)
setwd("/home/devdata/nyy/nyy_GFP/GFP/filter_vcfnew")
metadata<-read.delim(file ='sample_info.txt',header =T, quote = "",sep = "\t",stringsAsFactors=FALSE)
metadata$Individual[14:16]<-c("cloningW","cloningPF","cloningJR")

mydata<-read.delim(file ='merge_15sample_filter_add_4method_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
    mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$Individual %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
###########add sum 
data1<- mydata %>%
    dplyr::group_by(Individual) %>%
    dplyr::summarise(number = sum(number))

data1$base_change=rep('ALL',nrow(data1))
data1$Status<-case_when(
    data1$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    data1$Individual %in% c('NC1','NC2','NC3') ~ "NC",
    data1$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)

mydata<-rbind(data1,mydata)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata$base_change<-factor(mydata$base_change, levels=c('ALL',"A>G","C>T","A>C","A>T","C>G","C>A"))

mydata<-mydata[mydata$base_change=='ALL',]
mydata<-merge(mydata,metadata[,1:4],by='Individual')
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(mydata, aes(x=PCR, y=number, color=PCR)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of SNPs")+
    scale_y_continuous(breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = PCR,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum,hide.ns = F)
plotfile='./fig/PCR_barplot_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 9, height = 6)
plotfile='./fig/PCR_barplot_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 9, height = 6)

####################################################################
####################################################################
####################################################################
rm(list=ls())
library(reshape2)
library(ggplot2)
# options(scipen = 1)
setwd("/home/devdata/nyy/nyy_GFP/GFP/filter_vcfnew")
metadata<-read.delim(file ='sample_info.txt',header =T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-read.delim(file ='raw_uniq_number.txt',header =F, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('Individual','Number1')
mydata$Number<-mydata$Number1/100
mydata<-mydata[mydata$Individual !='07027',]
mydata<-merge(mydata,metadata,by='Individual',all = T)

breaks<-pretty(range(0,mydata$Number), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(data = mydata, mapping = aes(x = factor(Individual), y = Number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
    scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
    labs(x = "Individual",y = "Number of variants (x100)")+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20,face="bold"),)+
    guides(fill=FALSE)
plotfile='./fig/4method_raw_uniq_number_dp20_bar.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 8)

####PCR
library(ggplot2)
library(ggpubr)
# options(scipen = 9)
# breaks<-pretty(range(mydata$Number), 8)
# maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=PCR, y=Number, color=PCR)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of variants (x100)")+
    scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = PCR,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum,hide.ns = TRUE)
plotfile='./fig/4method_raw_uniq_number_dp20_PCR_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_raw_uniq_number_dp20_PCR_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

####SCNT
library(ggplot2)
library(ggpubr)
# options(scipen = 9)
mydata1<-mydata[is.na(mydata$SCNT)==FALSE,]
# breaks<-pretty(range(mydata1$Number), 8)
# maximum<- breaks[length(breaks)]
maximum1=40
p<-ggplot(mydata1, aes(x=SCNT, y=Number, color=SCNT)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of variants (x100)")+
    scale_y_continuous(limits = c(0,maximum1),breaks =breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = SCNT,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum1,hide.ns = TRUE)
plotfile='./fig/4method_raw_uniq_number_dp20_SCNT_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_raw_uniq_number_dp20_SCNT_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

####Development
library(ggplot2)
library(ggpubr)
# options(scipen = 9)

mydata1<-mydata[is.na(mydata$Dev)==FALSE,]
# breaks<-pretty(range(mydata1$Number), 8)
# maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=Dev, y=Number, color=Dev)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of variants (x100)")+
    scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = Dev,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum,hide.ns = TRUE)
plotfile='./fig/4method_raw_uniq_number_dp20_Development_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_raw_uniq_number_dp20_Development_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

####Status
library(ggplot2)
library(ggpubr)
# options(scipen = 9)
mydata1<-mydata[is.na(mydata$Status)==FALSE,]
# breaks<-pretty(range(mydata1$Number), 8)
# maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=Status, y=Number, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of variants (x100)")+
    scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum,hide.ns = TRUE)
plotfile='./fig/4method_raw_uniq_number_dp20_Status_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_raw_uniq_number_dp20_Status_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
