########
rm(list=ls())
library(ggplot2)
# install.packages('venn')
library(venn)
library(reshape2)
setwd("/home/devdata/nyy/nyy_GFP/GFP/filter_vcfnew")
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
    outname<-paste('./fig/',i,'_SNP_overlap_venn.pdf',sep='')
    pdf(outname, width = 8, height = 8)
    venn(x, snames = sname,ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8,ilcs=1.5)
    dev.off()
}


rm(list=ls())
library(reshape2)
library(ggplot2)
options(scipen = 9)
setwd("/home/devdata/nyy/nyy_GFP/GFP/filter_vcfnew")
metadata<-read.delim(file ='sample_info.txt',header =T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-read.delim(file ='4method_overlap_raw_count_dp20_filter.txt',header =F, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('Number','Individual')
mydata<-mydata[mydata$Individual !='07027',]
mydata<-merge(mydata,metadata,by='Individual',all = T)

breaks<-pretty(range(0,mydata$Number), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(data = mydata, mapping = aes(x = factor(Individual), y = Number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
    scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
    labs(x = "Individual",y = "Number of SNPs")+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20,face="bold"),)+
    guides(fill=FALSE)
plotfile='./fig/4method_overlap_raw_count_dp20_bar.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 8)

####PCR
library(ggplot2)
library(ggpubr)
options(scipen = 9)
breaks<-pretty(range(mydata$Number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=PCR, y=Number, color=PCR)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of SNPs")+
    scale_y_continuous(breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = PCR,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-400,hide.ns = TRUE)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_PCR_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_PCR_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

####SCNT
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata1<-mydata[is.na(mydata$SCNT)==FALSE,]
breaks<-pretty(range(mydata1$Number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=SCNT, y=Number, color=SCNT)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of SNPs")+
    scale_y_continuous(breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = SCNT,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-400,hide.ns = TRUE)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_SCNT_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_SCNT_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

####Development
library(ggplot2)
library(ggpubr)
options(scipen = 9)

mydata1<-mydata[is.na(mydata$Dev)==FALSE,]
breaks<-pretty(range(mydata1$Number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=Dev, y=Number, color=Dev)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of SNPs")+
    scale_y_continuous(breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = Dev,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-400,hide.ns = TRUE)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_Development_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_Development_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

####Status
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata1<-mydata[is.na(mydata$Status)==FALSE,]
breaks<-pretty(range(mydata1$Number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=Status, y=Number, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of SNPs")+
    scale_y_continuous(breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-400,hide.ns = TRUE)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_Status_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_raw_count_dp20_barplot_Status_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

############################################
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
options(scipen = 9)
mydata<-read.delim(file ='merge_15sample_filter_add_4method.txt',header = F, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('Individual','base_change','number')
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
mydata$Individual<-factor(mydata$Individual, levels=unique(mydata$Individual))

p<-ggplot(data = mydata, mapping = aes(x = factor(base_change), y = number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
    scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
    labs(x = "Base change",y = "Number of SNVs")+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_blank())
plotfile='./fig/merge_15sample_filter_add_4method_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)
plotfile='./fig/merge_15sample_filter_add_4method_count.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)

#####
mydata<-read.delim(file ='merge_15sample_filter_add_4method_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
mydata$Individual<-factor(mydata$Individual, levels=unique(mydata$Individual))                                                        
p<-ggplot(data = mydata, mapping = aes(x = factor(base_change), y = number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    labs(x = "Base change",y = "Number of SNVs")+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_blank())
plotfile='./fig/merge_15sample_filter_add_4method_count_merge.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)
plotfile='./fig/merge_15sample_filter_add_4method_count_merge.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)

#########
#########boxplot
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
options(scipen = 9)
mydata<-read.delim(file ='merge_15sample_filter_add_4method_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
    mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$Individual %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)

mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(mydata, aes(x=base_change, y=number, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    # stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-1,hide.ns = TRUE,method.args = list(alternative = "less"))
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-1,hide.ns = TRUE)

plotfile='./fig/base_change_barplot_PC_NC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 9, height = 6)
plotfile='./fig/base_change_barplot_PC_NC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 9, height = 6)

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
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(mydata, aes(x=base_change, y=number, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('ALL','A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-1,hide.ns = TRUE)
plotfile='./fig/base_change_barplot_PC_NC_v1_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/base_change_barplot_PC_NC_v1_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)

############only sum
mydata<-mydata[mydata$base_change=='ALL',]
options(scipen = 9)
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "",y = "Number of SNPs")+
    scale_y_continuous(breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-100,hide.ns = TRUE)
plotfile='./fig/4method_overlap_filter_sum_barplot_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_filter_sum_barplot_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

pcnumber<-c(data1$number[data1$Status=='PC'])
ncnumber<-c(data1$number[data1$Status=='NC'])
pcmean<-mean(pcnumber)
ncmean<-mean(ncnumber)
pcse<-sd(pcnumber)/sqrt(length(pcnumber))
ncse<-sd(ncnumber)/sqrt(length(ncnumber))


#########################
#########boxplot proportion
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_15sample_filter_add_4method_count_merge_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]

mydata$Status<-case_when(
    mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$Individual %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)

mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))

# mydata<-mydata[mydata$Individual !='PC3',]
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))

p<-ggplot(mydata, aes(x=base_change, y=number, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Base change",y = "Percent of SNVs (%)")+
    scale_y_continuous(limits = c(0,100),breaks = seq(0,100,10))+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=70, hide.ns = TRUE)
plotfile='./fig/barplot_PC_NC_proportion_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/barplot_PC_NC_proportion_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

# x<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='PC']
# y<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='NC']
# t.test(x,y)$p.value

######PC v2
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_15sample_filter_add_4method_count_merge_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
    mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$Individual %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)

mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
BEdata<-mydata[mydata$Status=='PC',]
BEdata$base_change<-factor(BEdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
# my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(BEdata, aes(x=base_change, y=number, color=base_change)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Base change",y = "Percent of SNVs (%)")+
    # scale_y_continuous(limits = c(0,100),breaks = seq(0,100,10))+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))+
	guides(color=FALSE)+
    stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5)

plotfile='./fig/base_change_barplot_onlyPC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/base_change_barplot_onlyPC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


######NC
CTdata<-mydata[mydata$Status=='NC',]
x<-CTdata$number[CTdata$base_change=="A>G" & CTdata$Status=='NC']
y<-CTdata$number[CTdata$base_change=="A>C" & CTdata$Status=='NC']
t.test(x,y)$p.value

CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
# my_comparisons <- list(c('A>G', 'C>T'),c('C>T', "A>C"),c('C>T', "A>T"),c('C>T', "C>G"),c('C>T',"C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(CTdata, aes(x=base_change, y=number, color=base_change)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Base change",y = "Percent of SNVs (%)")+
    scale_y_continuous(limits = c(0,100),breaks = seq(0,100,20))+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5)
plotfile='./fig/base_change_barplot_onlyNC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/base_change_barplot_onlyNC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

#############number
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_15sample_filter_add_4method_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
    mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$Individual %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))

BEdata<-mydata[mydata$Status=='PC',]
maximum<-max(BEdata$number)
maximum
breaks<-pretty(range(0,maximum*1.5), 8)
maximum<- breaks[length(breaks)]

BEdata$base_change<-factor(BEdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(BEdata, aes(x=base_change, y=number, color=base_change)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
    guides(color=FALSE)+
	stat_compare_means(comparisons=my_comparisons, method = "wilcox.test",label = "p.format",size=5)

plotfile='./fig/count_barplot_onlyPC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/count_barplot_onlyPC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

######NC
CTdata<-mydata[mydata$Status=='NC',]
CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))
range(CTdata$number) ## 13 343
breaks<-pretty(range(0,550), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(CTdata, aes(x=base_change, y=number, color=base_change)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=3)
plotfile='./fig/count_barplot_onlyNC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/count_barplot_onlyNC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

#########gene_region -------------------------------------------------------------
####all_change type
rm(list=ls())
library(tidyverse)
mydata<-read.delim(file ='SNP_gene_region_count.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
#mydata<-mydata[mydata$type %in% c('A>G','T>C'),]
mydata1 <- mydata %>%
    dplyr::select(sample, Region, count) %>%
    dplyr::group_by(sample,Region) %>% 
    dplyr::summarise(count = sum(count))
mydata1<-rbind(mydata1, data.frame(sample="NC2", Region = "UTR",count=0))

mydata2<- mydata1 %>%
    dplyr::group_by(sample) %>% 
    dplyr::mutate(proportion = count /sum(count))
mydata<-mydata2
mydata$Status<-case_when(
    mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$sample %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
maximum<-max(mydata$count)
maximum
breaks<-pretty(range(0,maximum), 6)
maximum<- breaks[length(breaks)]
mydata$Region<-factor(mydata$Region, levels=c("UTR","Exonic","Intronic","Downstream","Upstream","Intergenic"))

p<-ggplot(mydata, aes(x=Region, y=count, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Region",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-50,hide.ns = TRUE)
plotfile='./fig/SNP_gene_region_count_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/SNP_gene_region_count_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

breaks<-pretty(range(0,1), 5)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Region, y=proportion, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Region",y = "Percent of SNVs (%)")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=0.8,hide.ns = TRUE)

plotfile='./fig/SNP_gene_region_proportion_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/SNP_gene_region_proportion_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


####################################only c('A>G','T>C'),]
rm(list=ls())
library(tidyverse)
mydata<-read.delim(file ='SNP_gene_region_count.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
mydata<-mydata[mydata$type %in% c('A>G','T>C'),]
mydata1 <- mydata %>%
    dplyr::select(sample, Region, count) %>%
    dplyr::group_by(sample,Region) %>% 
    dplyr::summarise(count = sum(count))
tmp<-data.frame(sample=c("NC1","NC2","PC1","PC2","PC7"), Region = c("UTR","UTR","UTR","UTR","UTR"),count=c(0,0,0,0,0))
tmp1<-data.frame(sample=c("NC2","PC4"), Region = c("Exonic","Exonic"),count=c(0,0))
mydata1<-rbind(mydata1, tmp,tmp1)

mydata2<- mydata1 %>%
    dplyr::group_by(sample) %>% 
    dplyr::mutate(proportion = count /sum(count))
mydata<-mydata2
mydata$Status<-case_when(
    mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$sample %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata$Region<-factor(mydata$Region, levels=c("UTR","Exonic","Intronic","Downstream","Upstream","Intergenic"))
maximum<-max(mydata$count)
maximum
breaks<-pretty(range(0,maximum), 6)
maximum<- breaks[length(breaks)]

p<-ggplot(mydata, aes(x=Region, y=count, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Region",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=maximum-50,hide.ns = TRUE)
plotfile='./fig/SNP_gene_region_count_only_A2G_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/SNP_gene_region_count_only_A2G_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

breaks<-pretty(range(0,1), 5)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Region, y=proportion, color=Status)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Region",y = "Percent of SNVs (%)")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,label.y=0.8,hide.ns = TRUE)

plotfile='./fig/SNP_gene_region_proportion_only_A2G_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/SNP_gene_region_proportion_only_A2G_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

#########gene region onlyPC v1
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(tidyverse)l
mydata<-read.delim(file ='SNP_gene_region_count.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
#mydata<-mydata[mydata$type %in% c('A>G','T>C'),]
mydata1 <- mydata %>%
    dplyr::select(sample, Region, count) %>%
    dplyr::group_by(sample,Region) %>% 
    dplyr::summarise(count = sum(count))
mydata1<-rbind(mydata1, data.frame(sample="NC2", Region = "UTR",count=0))

mydata2<- mydata1 %>%
    dplyr::group_by(sample) %>% 
    dplyr::mutate(proportion = count /sum(count))
mydata<-mydata2
mydata$Status<-case_when(
    mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$sample %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata<-mydata[mydata$Status=='PC',]

breaks<-pretty(range(mydata$count), 8)
maximum<- breaks[length(breaks)]
mydata$Region<-factor(mydata$Region, levels=c("UTR","Exonic","Intronic","Downstream","Upstream","Intergenic"))

p<-ggbarplot(mydata, x = "Region", y = "count",size=1.5,color = "Region", ##fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Region",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    # scale_x_discrete(labels=c("UTR","CDS","Up_downstream","Intron","intergenic"))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    guides(color=FALSE)+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))
plotfile='offtargetRegion_barplot_onlyPC_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='offtargetRegion_barplot_onlyPC_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p<-ggplot(mydata, aes(x=Region, y=count, color=Region)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Region",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    guides(color=FALSE)
plotfile='./fig/offtargetRegion_barplot_onlyPC_v2.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/offtargetRegion_barplot_onlyPC_v2.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

#########################gene region onlyPC and  only c('A>G','T>C'),]
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(tidyverse)
mydata<-read.delim(file ='SNP_gene_region_count.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
mydata<-mydata[mydata$type %in% c('A>G','T>C'),]
mydata1 <- mydata %>%
    dplyr::select(sample, Region, count) %>%
    dplyr::group_by(sample,Region) %>% 
    dplyr::summarise(count = sum(count))
tmp<-data.frame(sample=c("NC1","NC2","PC1","PC2","PC7"), Region = c("UTR","UTR","UTR","UTR","UTR"),count=c(0,0,0,0,0))
tmp1<-data.frame(sample=c("NC2","PC4"), Region = c("Exonic","Exonic"),count=c(0,0))
mydata1<-rbind(mydata1, tmp,tmp1)

mydata2<- mydata1 %>%
    dplyr::group_by(sample) %>% 
    dplyr::mutate(proportion = count /sum(count))
mydata<-mydata2
mydata$Status<-case_when(
    mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "PC",
    mydata$sample %in% c('NC1','NC2','NC3') ~ "NC",
    mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "PC"
)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata<-mydata[mydata$Status=='PC',]

breaks<-pretty(range(mydata$count), 8)
maximum<- breaks[length(breaks)]
mydata$Region<-factor(mydata$Region, levels=c("UTR","Exonic","Intronic","Downstream","Upstream","Intergenic"))

p<-ggbarplot(mydata, x = "Region", y = "count",size=1.5,color = "Region", ##fill="Status",
             add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Region",y = "Number of SNVs")+
    scale_y_continuous(breaks = breaks)+
    # scale_x_discrete(labels=c("UTR","CDS","Up_downstream","Intron","intergenic"))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    guides(color=FALSE)+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))
plotfile='./fig/offtargetRegion_barplot_onlyPC_onlyA2G_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/offtargetRegion_barplot_onlyPC_onlyA2G_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p<-ggplot(mydata, aes(x=Region, y=count, color=Region)) +
    geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
    geom_jitter(position=position_dodge(0.8), size=0.8)+
    labs(x = "Region",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    guides(color=FALSE)
plotfile='./fig/offtargetRegion_barplot_onlyPC_onlyA2G_v2.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/offtargetRegion_barplot_onlyPC_onlyA2G_v2.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


############ggseqlogo
rm(list=ls())
# Load the required packages
library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(stringi)
library(patchwork)
mydata<-read.delim(file ='all_edit_site_seqlog.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq','changetype','genetype','sample')
# mydata<-mydata[mydata$changetype=='A>G',]
# mydata<-mydata[mydata$changetype=='T>C',]
TC_seq<-mydata$seq[mydata$changetype=='T>C']
AG_seq<-mydata$seq[mydata$changetype=='A>G']

rev_seq <- c(stri_reverse(TC_seq),AG_seq)
p1<-ggseqlogo(rev_seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 

all_seq<-c(TC_seq,AG_seq)
p2<-ggseqlogo(all_seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
p1+p2
plotfile='./fig/all_edit_ggseqlogo.pdf'
ggsave(plotfile, plot=p2, dpi = 600,width = 8, height = 6)
plotfile='./fig/all_edit_ggseqlogo.png'
ggsave(plotfile, plot=p2, dpi = 600,width = 8, height = 6)


p1mydata<-mydata[mydata$sample=='PC1',]
p<-ggseqlogo(p1mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC1.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC1.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p2mydata<-mydata[mydata$sample=='PC2',]
p<-ggseqlogo(p2mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC2.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC2.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p3mydata<-mydata[mydata$sample=='PC3',]
p<-ggseqlogo(p3mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC3.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC3.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p4mydata<-mydata[mydata$sample=='PC4',]
p<-ggseqlogo(p4mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC4.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC4.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p5mydata<-mydata[mydata$sample=='PC5',]
p<-ggseqlogo(p6mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC5.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC5.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p6mydata<-mydata[mydata$sample=='PC6',]
p<-ggseqlogo(p6mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC6.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC6.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


p7mydata<-mydata[mydata$sample=='PC7',]
p<-ggseqlogo(p7mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC7.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC7.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

p8mydata<-mydata[mydata$sample=='PC8',]
p<-ggseqlogo(p8mydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/PC8.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC8.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

wmydata<-mydata[mydata$sample=='cloningW',]
p<-ggseqlogo(wmydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/cloningW.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/cloningW.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

PFmydata<-mydata[mydata$sample=='cloningPF',]
p<-ggseqlogo(PFmydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/cloningPF.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/cloningPF.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

JRmydata<-mydata[mydata$sample=='cloningJR',]
p<-ggseqlogo(JRmydata$seq,method='p')+ theme(axis.text=element_text(size=15,colour='black'),axis.title=element_text(size=16)) 
plotfile='./fig/cloningJR.genome_ggseqlogo.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/cloningJR.genome_ggseqlogo.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

