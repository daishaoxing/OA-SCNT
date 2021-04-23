rm(list=ls())
library(reshape2)
library(ggplot2)
options(scipen = 9)

mydata<-read.delim(file ='GATK_10sample_offtarget_count.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('Individual','base_change','number')
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(data = mydata, mapping = aes(x = factor(base_change), y = number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
    scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
    labs(x = "Base change",y = "Number of SNVs")+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_text(size=20))
plotfile='GATK_10sample_offtarget_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 5)

#####
mydata<-read.delim(file ='GATK_10sample_offtarget_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('Individual','base_change','number')
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
                                                        
p<-ggplot(data = mydata, mapping = aes(x = factor(base_change), y = number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    labs(x = "Base change",y = "Number of SNVs")+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_text(size=20))
plotfile='GATK_10sample_offtarget_count_merge.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 5)

#########
#########boxplot
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='GATK_10sample_offtarget_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
Status<-c()
a<-rep('NC',nrow(mydata[mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5'),]))
b<-rep('PC',nrow(mydata[mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5'),]))
mydata$Status<-c(a,b)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]+100
p<-ggbarplot(mydata, x = "base_change", y = "number",size=1.5,color = "Status",# fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title =element_blank())+
    stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=maximum-1)
plotfile='GATKbase_change_barplot_PC_NCofftarget.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='GATKbase_change_barplot_PC_NCofftarget.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

# x<-mydata$number[mydata$base_change=="G>A" & mydata$Status=='PC']
# y<-mydata$number[mydata$base_change=="G>A" & mydata$Status=='NC']
# t.test(x,y)$p.value
# x<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='PC']
# y<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='NC']
# t.test(x,y)$p.value

# x<-mydata$number[mydata$base_change=="G>A" & mydata$Status=='NC']
# y<-mydata$number[mydata$base_change=="G>A" & mydata$Status=='NC']
# t.test(x,y)$p.value
# x<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='PC']
# y<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='PC']
# t.test(x,y)$p.value

###########add sum 
number<-c(1:10)
number[1]<-sum(mydata$number[mydata$Individual=='NC1' & mydata$Status=='NC'])
number[2]<-sum(mydata$number[mydata$Individual=='NC2' & mydata$Status=='NC'])
number[3]<-sum(mydata$number[mydata$Individual=='NC3' & mydata$Status=='NC'])
number[4]<-sum(mydata$number[mydata$Individual=='NC4' & mydata$Status=='NC'])
number[5]<-sum(mydata$number[mydata$Individual=='NC5' & mydata$Status=='NC'])
number[6]<-sum(mydata$number[mydata$Individual=='PC1' & mydata$Status=='PC'])
number[7]<-sum(mydata$number[mydata$Individual=='PC2' & mydata$Status=='PC'])
number[8]<-sum(mydata$number[mydata$Individual=='PC3' & mydata$Status=='PC'])
number[9]<-sum(mydata$number[mydata$Individual=='PC4' & mydata$Status=='PC'])
number[10]<-sum(mydata$number[mydata$Individual=='PC5' & mydata$Status=='PC'])

data1<-data.frame(Individual=c('NC1','NC2','NC3','NC4','NC5','PC1','PC2','PC3','PC4','PC5'),base_change=rep('ALL',10),number=number,Status=c(rep('NC',5),rep('PC',5)))

mydata<-rbind(data1,mydata)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata$base_change<-factor(mydata$base_change, levels=c('ALL',"A>G","C>T","A>C","A>T","C>G","C>A"))
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
p<-ggbarplot(mydata, x = "base_change", y = "number",size=1.5,color = "Status", ##fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('ALL','A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=maximum-1)
plotfile='GATKbase_change_barplot_PC_NCofftargetv1a.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='GATKbase_change_barplot_PC_NCofftargetv1a.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


#########################
#########boxplot proportion
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='GATK_10sample_offtarget_count_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
Status<-c()
a<-rep('NC',nrow(mydata[mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5'),]))
b<-rep('PC',nrow(mydata[mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5'),]))
mydata$Status<-c(a,b)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
p<-ggbarplot(mydata, x = "base_change", y = "number",size=1.5,color = "Status", ##fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Base change",y = "Percent of SNVs (%)")+
    scale_y_continuous(limits = c(0,60),breaks = seq(0,60,10))+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
    stat_compare_means(aes(group = Status),method = "t.test",label = "p.signif",size=8,label.y=55)
plotfile='GATKbarplot_PC_NC_proportionofftargeta.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='GATKbarplot_PC_NC_proportionofftargeta.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

x<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='PC']
y<-mydata$number[mydata$base_change=="A>G" & mydata$Status=='NC']
t.test(x,y)$p.value
x<-mydata$number[mydata$base_change=="C>T" & mydata$Status=='PC']
y<-mydata$number[mydata$base_change=="C>T" & mydata$Status=='NC']
t.test(x,y)$p.value


######PC v2
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='GATK_10sample_offtarget_count_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
Status<-c()
a<-rep('NC',nrow(mydata[mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5'),]))
b<-rep('PC',nrow(mydata[mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5'),]))
mydata$Status<-c(a,b)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
BEdata<-mydata[mydata$Status=='PC',]
BEdata$base_change<-factor(BEdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
# my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))
p<-ggbarplot(BEdata, x = "base_change", y = "number",size=1.5,color = "base_change", ##fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Base change",y = "Percent of SNVs (%)")+
    scale_y_continuous(limits = c(0,80),breaks = seq(0,80,10))+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=3)
plotfile='GATK_base_change_barplot_onlyPCofftarget.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='GATK_base_change_barplot_onlyPCofftarget.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


######NC

CTdata<-mydata[mydata$Status=='NC',]
x<-CTdata$number[CTdata$base_change=="A>G" & CTdata$Status=='NC']
y<-CTdata$number[CTdata$base_change=="A>G" & CTdata$Status=='NC']
t.test(x,y)$p.value

CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggbarplot(CTdata, x = "base_change", y = "number",size=1.5,color = "base_change", ##fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Base change",y = "Percent of SNVs (%)")+
    scale_y_continuous(limits = c(-1,60),breaks = seq(0,60,10))+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=3)
plotfile='GATK_base_change_barplot_onlyNCofftargeta.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='GATK_base_change_barplot_onlyNCofftargeta.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


#############number
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='GATK_10sample_offtarget_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
Status<-c()
a<-rep('NC',nrow(mydata[mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5'),]))
b<-rep('PC',nrow(mydata[mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5'),]))
mydata$Status<-c(a,b)
mydata$Status<-factor(mydata$Status, levels=c('PC', 'NC'))
BEdata<-mydata[mydata$Status=='PC',]
range(BEdata$number) ##132 1797
breaks<-pretty(range(0,3000), 8)
maximum<- breaks[length(breaks)]
BEdata$base_change<-factor(BEdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
# my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))
p<-ggbarplot(BEdata, x = "base_change", y = "number",size=1.5,color = "base_change", ##fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=3)
plotfile='GATK_count_barplot_onlyPCofftargeta.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='GATK_count_barplot_onlyPCofftargeta.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

######NC

CTdata<-mydata[mydata$Status=='NC',]

# x<-CTdata$number[CTdata$base_change=="A>G" & CTdata$Status=='NC']
# y<-CTdata$number[CTdata$base_change=="A>G" & CTdata$Status=='NC']
# t.test(x,y)$p.value

CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))
range(CTdata$number) ## 7 108
breaks<-pretty(range(0,180), 8)
maximum<- breaks[length(breaks)]
p<-ggbarplot(CTdata, x = "base_change", y = "number",size=1.5,color = "base_change", ##fill="Status",
    add = c("jitter",'mean_se'),merge=TRUE)+
    labs(x = "Base change",y = "Number of SNVs")+
    scale_y_continuous(limits = c(-1,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
    guides(color=FALSE)+
    stat_compare_means(comparisons=my_comparisons,method = "t.test",label = "p.signif",size=3)
plotfile='GATK_count_barplot_onlyNCofftargeta.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='GATK_count_barplot_onlyNCofftargeta.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

#########jitterplot
setwd("/home/devdata/nyy/nyy_GFP/RNA_offtarget")
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_GATK_10sample_Jitter_v1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
##"Individual" "percent"    "Status"     "chr"
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
p<-ggplot(data = mydata, mapping = aes(x = Individual, y = percent)) + geom_jitter(colour='red',width = 0.4,size = 1.5)+
scale_y_continuous(limits = c(0,101),breaks =breaks)+
labs(x = '',y = "RNA\nA-to-G editing (%)")+
guides(color=FALSE)+
theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))

plotfile='./fig/merge_GATK_10sample_Jitter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/merge_GATK_10sample_Jitter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)


############################V2#################
mydata<-read.delim(file ='merge_GATK_10sample_Jitter_v1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
##"Individual" "percent"    "Status"     "chr"
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
p<-ggplot(data = mydata, mapping = aes(x = Individual, y = percent,color=Status)) + geom_jitter(width = 0.4,size = 0.001)+
  scale_y_continuous(limits = c(0,101),breaks =breaks)+
  labs(x = '',y = "RNA\nA-to-G editing (%)")+
  guides(color=FALSE)+
  scale_color_manual(values =c('#0571b0','#ca0020'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))

plotfile='./fig/merge_GATK_10sample_Jitter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/merge_GATK_10sample_Jitter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)

chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y")

####chr
mydata<-read.delim(file ='merge_GATK_10sample_Jitter_v1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[mydata$Individual=='PC1' & mydata$chr %in% chrs,]
breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
mydata$chr<-factor(mydata$chr, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))
p<-ggplot(data = mydata, mapping = aes(x = chr, y = percent,colour=chr)) + geom_jitter(width = 0.4,size = 1.5)+
scale_y_continuous(limits = c(0,101),breaks =breaks)+
labs(x = '',y = "RNA\nA-to-G editing (%)")+
scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))+
theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
guides(color=FALSE)+
theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))

plotfile='./fig/GATK_PC1sample_Jitter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/GATK_PC1sample_Jitter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)

####chr

mydata<-read.delim(file ='merge_GATK_10sample_Jitter_v1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[mydata$Individual=='PC2' & mydata$chr %in% chrs,]

breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
mydata$chr<-factor(mydata$chr, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))
p<-ggplot(data = mydata, mapping = aes(x = chr, y = percent,colour=chr)) + geom_jitter(width = 0.4,size = 1.5)+
scale_y_continuous(limits = c(0,101),breaks =breaks)+
labs(x = '',y = "RNA\nA-to-G editing (%)")+
scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))+
theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
guides(color=FALSE)+
theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))

plotfile='./fig/GATK_PC2sample_Jitter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/GATK_PC2sample_Jitter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)

####chr
mydata<-read.delim(file ='merge_GATK_10sample_Jitter_v1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

mydata<-mydata[mydata$Individual=='PC3' & mydata$chr %in% chrs,]
breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
mydata$chr<-factor(mydata$chr, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))
p<-ggplot(data = mydata, mapping = aes(x = chr, y = percent,colour=chr)) + geom_jitter(width = 0.4,size = 1.5)+
scale_y_continuous(limits = c(0,101),breaks =breaks)+
labs(x = '',y = "RNA\nA-to-G editing (%)")+
scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))+
theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
guides(color=FALSE)+
theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))
plotfile='./fig/GATK_PC3sample_Jitter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/GATK_PC3sample_Jitter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)

####chr
mydata<-read.delim(file ='merge_GATK_10sample_Jitter_v1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

mydata<-mydata[mydata$Individual=='PC4' & mydata$chr %in% chrs,]
breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
mydata$chr<-factor(mydata$chr, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))
p<-ggplot(data = mydata, mapping = aes(x = chr, y = percent,colour=chr)) + geom_jitter(width = 0.4,size = 1.5)+
scale_y_continuous(limits = c(0,101),breaks =breaks)+
labs(x = '',y = "RNA\nA-to-G editing (%)")+
scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))+
theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
guides(color=FALSE)+
theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))
plotfile='./fig/GATK_PC4sample_Jitter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/GATK_PC4sample_Jitter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)

####chr
mydata<-read.delim(file ='merge_GATK_10sample_Jitter_v1.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

mydata<-mydata[mydata$Individual=='PC5' & mydata$chr %in% chrs,]
breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
mydata$chr<-factor(mydata$chr, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))
p<-ggplot(data = mydata, mapping = aes(x = chr, y = percent,colour=chr)) + geom_jitter(width = 0.4,size = 1.5)+
scale_y_continuous(limits = c(0,101),breaks =breaks)+
labs(x = '',y = "RNA\nA-to-G editing (%)")+
scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))+
theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
guides(color=FALSE)+
theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))
plotfile='./fig/GATK_PC5sample_Jitter_v1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/GATK_PC5sample_Jitter_v1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)


rm(list=ls())
# ggseqlogo######Load the required packages
require(ggplot2)
require(ggseqlogo)
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='all_RNAedit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='all_RNAedit_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='all_RNAedit_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC1.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC1.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC1.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC2.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC2.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC2.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC3.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC3.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC3.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC4.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC4.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC4.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC5.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC5.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC5.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC6.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC6.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC6.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC7.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC7.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC7.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC8.RNA_edit_siteseqAG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC8.RNA_ggseqlogoAG.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC8.RNA_ggseqlogoAG.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
###################CT
mydata<-read.delim(file ='all_RNAedit_siteseqCT.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='all_RNAedit_ggseqlogoCT.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='all_RNAedit_ggseqlogoCT.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC1.RNA_edit_siteseqCT.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC1.RNA_ggseqlogoCT.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC1.RNA_ggseqlogoCT.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC2.RNA_edit_siteseqCT.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC2.RNA_ggseqlogoCT.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC2.RNA_ggseqlogoCT.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC3.RNA_edit_siteseqCT.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC3.RNA_ggseqlogoCT.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC3.RNA_ggseqlogoCT.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC4.RNA_edit_siteseqCT.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC4.RNA_ggseqlogoCT.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC4.RNA_ggseqlogoCT.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
mydata<-read.delim(file ='PC5.RNA_edit_siteseqCT.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq')
p<-ggseqlogo(mydata$seq,method='p')
plotfile='PC5.RNA_ggseqlogoCT.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC5.RNA_ggseqlogoCT.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)


################gene expression
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(reshape2)
options(scipen = 9)
mygene<-read.table(file ='A2Ggene.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mygene)<-c('gene','sample')
data_exp<-read.table(file ='./3_RNA_analysis_output/TPM.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

colnames(data_exp)<-c('gene',sapply(strsplit(colnames(data_exp),".",fixed = TRUE), `[`, 2)[-1])
rownames(data_exp)<-data_exp$gene
data_exp<-log2(data_exp[,-1]+1)

########PC1
set.seed(1234)
offtarget<-mygene$gene[mygene$sample=='PC1']
genes<-setdiff(rownames(data_exp),offtarget)
random<-sample(genes,length(offtarget))
PC1data<-data.frame(PC1=data_exp$PC1[rownames(data_exp) %in% offtarget],
                   random_genes=data_exp$PC1[rownames(data_exp) %in% random],
				   NC1=data_exp$NC1[rownames(data_exp) %in% offtarget],
                   NC2=data_exp$NC2[rownames(data_exp) %in% offtarget],
				   NC3=data_exp$NC3[rownames(data_exp) %in% offtarget],
				   NC4=data_exp$NC4[rownames(data_exp) %in% offtarget],
				   NC5=data_exp$NC5[rownames(data_exp) %in% offtarget])
PC1data <- melt(PC1data,  variable.name = "sample", value.name = "TPM")
breaks<-pretty(range(PC1data$TPM), 8)
maximum<- breaks[length(breaks)]+2
p<-ggplot(PC1data,aes(x=sample, y=TPM, fill=sample)) +
    geom_boxplot(outlier.colour = NA) +
    # geom_jitter(size=0.4,position=position_jitter(0.1)) +
    # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
    labs(x = "",y = "log2(TPM+1)")+
    scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('PC1 offtarget\ngenes', 'PC1 Random\ngenes', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    guides(fill=FALSE)+
	theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title=element_text(size=20))+
    stat_compare_means(ref.group = "PC1",method = "t.test",label = "p.signif",size=8,label.y=maximum-1)
plotfile='PC1_RNA_offtarget_exp.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC1_RNA_offtarget_exp.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

########PC2
set.seed(1234)
offtarget<-mygene$gene[mygene$sample=='PC2']
genes<-setdiff(rownames(data_exp),offtarget)
random<-sample(genes,length(offtarget))
PC2data<-data.frame(PC2=data_exp$PC2[rownames(data_exp) %in% offtarget],
                   random_genes=data_exp$PC2[rownames(data_exp) %in% random],
				   NC1=data_exp$NC1[rownames(data_exp) %in% offtarget],
                   NC2=data_exp$NC2[rownames(data_exp) %in% offtarget],
				   NC3=data_exp$NC3[rownames(data_exp) %in% offtarget],
				   NC4=data_exp$NC4[rownames(data_exp) %in% offtarget],
				   NC5=data_exp$NC5[rownames(data_exp) %in% offtarget])
PC2data <- melt(PC2data,  variable.name = "sample", value.name = "TPM")
breaks<-pretty(range(PC2data$TPM), 8)
maximum<- breaks[length(breaks)]+2
p<-ggplot(PC2data,aes(x=sample, y=TPM, fill=sample)) +
    geom_boxplot(outlier.colour = NA) +
    # geom_jitter(size=0.4,position=position_jitter(0.1)) +
    # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
    labs(x = "",y = "log2(TPM+1)")+
    scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('PC2 offtarget\ngenes', 'PC2 Random\ngenes', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    guides(fill=FALSE)+
	theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title=element_text(size=20))+
    stat_compare_means(ref.group = "PC2",method = "t.test",label = "p.signif",size=8,label.y=maximum-1)
plotfile='PC2_RNA_offtarget_exp.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC2_RNA_offtarget_exp.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
########PC3
set.seed(1234)
offtarget<-mygene$gene[mygene$sample=='PC3']
genes<-setdiff(rownames(data_exp),offtarget)
random<-sample(genes,length(offtarget))
PC3data<-data.frame(PC3=data_exp$PC3[rownames(data_exp) %in% offtarget],
                   random_genes=data_exp$PC3[rownames(data_exp) %in% random],
				   NC1=data_exp$NC1[rownames(data_exp) %in% offtarget],
                   NC2=data_exp$NC2[rownames(data_exp) %in% offtarget],
				   NC3=data_exp$NC3[rownames(data_exp) %in% offtarget],
				   NC4=data_exp$NC4[rownames(data_exp) %in% offtarget],
				   NC5=data_exp$NC5[rownames(data_exp) %in% offtarget])
PC3data <- melt(PC3data,  variable.name = "sample", value.name = "TPM")
breaks<-pretty(range(PC3data$TPM), 8)
maximum<- breaks[length(breaks)]+2
p<-ggplot(PC3data,aes(x=sample, y=TPM, fill=sample)) +
    geom_boxplot(outlier.colour = NA) +
    # geom_jitter(size=0.4,position=position_jitter(0.1)) +
    # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
    labs(x = "",y = "log2(TPM+1)")+
    scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('PC3 offtarget\ngenes', 'PC3 Random\ngenes', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    guides(fill=FALSE)+
	theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title=element_text(size=20))+
    stat_compare_means(ref.group = "PC3",method = "t.test",label = "p.signif",size=8,label.y=maximum-1)
plotfile='PC3_RNA_offtarget_exp.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC3_RNA_offtarget_exp.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
########PC4
set.seed(1234)
offtarget<-mygene$gene[mygene$sample=='PC4']
genes<-setdiff(rownames(data_exp),offtarget)
random<-sample(genes,length(offtarget))
PC4data<-data.frame(PC4=data_exp$PC4[rownames(data_exp) %in% offtarget],
                   random_genes=data_exp$PC4[rownames(data_exp) %in% random],
				   NC1=data_exp$NC1[rownames(data_exp) %in% offtarget],
                   NC2=data_exp$NC2[rownames(data_exp) %in% offtarget],
				   NC3=data_exp$NC3[rownames(data_exp) %in% offtarget],
				   NC4=data_exp$NC4[rownames(data_exp) %in% offtarget],
				   NC5=data_exp$NC5[rownames(data_exp) %in% offtarget])
PC4data <- melt(PC4data,  variable.name = "sample", value.name = "TPM")
breaks<-pretty(range(PC4data$TPM), 8)
maximum<- breaks[length(breaks)]+2
p<-ggplot(PC4data,aes(x=sample, y=TPM, fill=sample)) +
    geom_boxplot(outlier.colour = NA) +
    # geom_jitter(size=0.4,position=position_jitter(0.1)) +
    # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
    labs(x = "",y = "log2(TPM+1)")+
    scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('PC4 offtarget\ngenes', 'PC4 Random\ngenes', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    guides(fill=FALSE)+
	theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title=element_text(size=20))+
    stat_compare_means(ref.group = "PC4",method = "t.test",label = "p.signif",size=8,label.y=maximum-1)
plotfile='PC4_RNA_offtarget_exp.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC4_RNA_offtarget_exp.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
########PC5
set.seed(1234)
offtarget<-mygene$gene[mygene$sample=='PC5']
genes<-setdiff(rownames(data_exp),offtarget)
random<-sample(genes,length(offtarget))
PC5data<-data.frame(PC5=data_exp$PC5[rownames(data_exp) %in% offtarget],
                   random_genes=data_exp$PC5[rownames(data_exp) %in% random],
				   NC1=data_exp$NC1[rownames(data_exp) %in% offtarget],
                   NC2=data_exp$NC2[rownames(data_exp) %in% offtarget],
				   NC3=data_exp$NC3[rownames(data_exp) %in% offtarget],
				   NC4=data_exp$NC4[rownames(data_exp) %in% offtarget],
				   NC5=data_exp$NC5[rownames(data_exp) %in% offtarget])
PC5data <- melt(PC5data,  variable.name = "sample", value.name = "TPM")
breaks<-pretty(range(PC5data$TPM), 8)
maximum<- breaks[length(breaks)]+2
p<-ggplot(PC5data,aes(x=sample, y=TPM, fill=sample)) +
    geom_boxplot(outlier.colour = NA) +
    # geom_jitter(size=0.4,position=position_jitter(0.1)) +
    # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
    labs(x = "",y = "log2(TPM+1)")+
    scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
    scale_x_discrete(labels=c('PC5 offtarget\ngenes', 'PC5 Random\ngenes', 'NC1', 'NC2', 'NC3', 'NC4', 'NC5'))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    guides(fill=FALSE)+
	theme(axis.text=element_text(size=16),axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title=element_text(size=20))+
    stat_compare_means(ref.group = "PC5",method = "t.test",label = "p.signif",size=8,label.y=maximum-1)
plotfile='PC5_RNA_offtarget_exp.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='PC5_RNA_offtarget_exp.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

##############################venn
rm(list=ls())
library(ggplot2)
library(venn)
library(reshape2)
options(scipen = 9)
mygene<-read.table(file ='A2Ggene.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mygene)<-c('gene','sample')
pc1 = mygene$gene[mygene$sample == 'PC1'] #??ͬ sample ?? gene ?ǹ??ൽһ?? list?У?????
pc2 = mygene$gene[mygene$sample == 'PC2']
pc3 = mygene$gene[mygene$sample == 'PC3']
pc4 = mygene$gene[mygene$sample == 'PC4']
pc5 = mygene$gene[mygene$sample == 'PC5']
c(length(pc1),length(pc2),length(pc3),length(pc4),length(pc5))
# 950 970 960 808 999
x = list(pc1,pc2,pc3,pc4,pc5)   #??????????Ҳ???????? factor ????һ?? list ??
# venn(x, snames = "PC1,PC2,PC3,PC4,PC5", zcolor = "style", cexil = 1, cexsn = 0.8)
tiff("5PC_gene_venn.tif", res = 600,width = 6, height = 6, units = 'in')
venn(x, snames = "PC1(950 genes),PC2(970 genes),PC3(960 genes),PC4(808 genes),PC5(999 genes)", ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8)
dev.off()

####KEGG enrichment-------------------------------------------------------------------
rm(list=ls())
library( "clusterProfiler")
library("org.Mmu.eg.db")
setwd('/home/devdata/nyy/nyy_GFP/GFP/integrated')
columns(org.Mmu.eg.db)
 # [1] "ACCNUM"       "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
 # [6] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"          
# [11] "GOALL"        "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PMID"        
# [16] "REFSEQ"       "SYMBOL"       "UNIPROT"     
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c("ENSEMBL",'ENTREZID'))
mygene<-read.table(file ='A2Ggene.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mygene)<-c('ENSEMBL','sample')
merger_gene<-merge(mygene,geneinfo,by='ENSEMBL',all.x = TRUE)
merger_gene<-na.omit(merger_gene)
write.table(merger_gene, file = "offtarget_merger_gene.tab", quote = FALSE,sep="\t",row.names = FALSE)   #208
####
mygene<-read.table(file ='offtarget_merger_gene.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ekegg <- enrichKEGG(mygene$ENTREZID,organism="mcc",pvalueCutoff=0.05)
ekegg<-as.data.frame(ekegg@result)
write.table(ekegg, file = "sig_kegg.tab", quote = FALSE,sep="\t",row.names = FALSE)   #208

# 8.bubble-plot KEGG----------------------------------------------------------------------

rm(list=ls())
library(ggplot2)
res <- c('sig_kegg')

for (i in res) {
  pathway = pathway<-read.delim(file = paste(i,'.tab',sep=''),header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  p <- ggplot(pathway[1:20,],aes(x=p.adjust,y=reorder(Description,-p.adjust),size=Count,color=p.adjust,ylab=''))+
    geom_point()+scale_colour_gradient(low="blue",high="red")+ 
    scale_size_area(name="genecounts")+theme_bw()+
    scale_size_continuous(range=c(4,10))+
    labs(y='',x='p.adjust')+
    geom_vline(xintercept = 0.05,linetype =2,colour = 'black')+
    theme(axis.text=element_text(size=15,face = "bold", color="gray50"),
          axis.title=element_text(size=20),
          legend.text=element_text(size=16),legend.title = element_text(size=20))
  plotfile = paste(i,'_pathway.png')
  ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 8)
}

###########GO
# BP
ego_BP <- enrichGO(gene =mygene$ENTREZID,
                   OrgDb=org.Mmu.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP_result<-as.data.frame(ego_BP@result)#1895
write.table(ego_BP_result, file = "GO_BP_sig.tab", quote = FALSE,sep="\t",row.names = FALSE)
#MF
ego_MF <- enrichGO(gene =mygene$ENTREZID,
                   OrgDb=org.Mmu.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05)
ego_MF_result<-as.data.frame(ego_MF@result)#464
write.table(ego_MF_result, file = "GO_MF_sig.tab", quote = FALSE,sep="\t",row.names = FALSE)
#CC
ego_CC <- enrichGO(gene =mygene$ENTREZID,
                     OrgDb=org.Mmu.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05)
ego_CC_result<-as.data.frame(ego_CC@result)#338
write.table(ego_CC_result, file = "GO_CC_sig.tab", quote = FALSE,sep="\t",row.names = FALSE)


# 10.GO histone-plot ------------------------------------------------------
rm(list=ls())
library(ggplot2)

BP_res <- read.delim(file = 'GO_BP_sig.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
BP_res$group = rep('BP',nrow(BP_res));sum(BP_res$pvalue < 0.05)#148
BP_res = BP_res[1:10,]
MF_res <- read.delim(file = 'GO_MF_sig.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
MF_res$group = rep('MF',nrow(MF_res));sum(MF_res$pvalue < 0.05)#54
MF_res = MF_res[1:10,]
CC_res <- read.delim(file = 'GO_CC_sig.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
CC_res$group = rep('CC',nrow(CC_res));sum(CC_res$pvalue < 0.05)#54
CC_res = CC_res[1:10,]

top10_go = rbind(BP_res,MF_res,CC_res)
top10_go[4,'Description'] <- 'negative regulation of protein modification'
top10_go$group = factor(top10_go$group, levels=c('BP','MF','CC'))
top10_go$number = factor(rev(1:nrow(top10_go)))

p <- ggplot(data=top10_go, aes(x=number, y=Count, fill=group)) +
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() +
  theme_bw() + 
  scale_x_discrete(labels=rev(top10_go$Description)) +
  xlab("") +
  theme(axis.text=element_text(size=12,face = "bold", color="gray50"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),legend.title = element_blank()) +
  scale_fill_discrete(name = 'Group')
# + labs(title = "Top10 Enriched GO Terms")

p
plotfile='offtarget_gene_GOplot.png'
ggsave(plotfile, plot=p, dpi = 600,width = 14, height = 10)

