setwd('/home/devdata/nyy/nyy_GFP/RNA_offtarget/2021v3')
rm(list=ls())
library(reshape2)
library(ggplot2)
options(scipen = 9)
#####merge_GATK_10sample_count.txt
mydata<-read.delim(file ='merge_GATK_10sample_count.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('Individual','base_change','number')
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(data = mydata, mapping = aes(x = factor(base_change), y = number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
  scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
  labs(x = "Base change",y = "SNVs")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_text(size=20))
p
plotfile='./fig/merge_GATK_10sample_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 5)

#####merge_GATK_10sample_count_merge.txt
mydata<-read.delim(file ='merge_GATK_10sample_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
# colnames(mydata)<-c('Individual','base_change','number')
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(data = mydata, mapping = aes(x = factor(base_change), y = number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  labs(x = "Base change",y = "SNVs")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_text(size=20))
p
plotfile='./fig/merge_GATK_10sample_count_merge.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 5)

# change count heatmap ----------------------------------------------------
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
options(scipen = 9)
mydata<-read.delim(file ='merge_GATK_10sample_count.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('Individual','base_change','number')
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5') ~ "SCNT-ABE"
)
mydata<- mydata %>%
  dplyr::group_by(base_change,Status) %>%
  dplyr::mutate(mean_count = round(mean(number),0))

mydata$From <- str_split_fixed(mydata$base_change,"[>]",n=2)[,1]
mydata$To <- str_split_fixed(mydata$base_change,"[>]",n=2)[,2]
mydata<-unique(mydata[,4:7])

colnames(mydata)
tmpdf<-data.frame("Status"=c(rep('SCNT',4),rep('SCNT-ABE',4)),
                  "mean_count"=c(rep(0,8)),
                  "From"=c('A','T','G','C','A','T','G','C'),
                  "To"=c('A','T','G','C','A','T','G','C'))
mydata<-rbind(as.data.frame(mydata),tmpdf)
mydata$Status<-factor(mydata$Status, levels=c('SCNT', 'SCNT-ABE'))

p<-ggplot(mydata, aes(x = To, y=From)) +
  geom_tile(aes(fill = mean_count),colour = "black") +
  geom_text(aes(label = mean_count),size=8,colour='black') +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_continuous(low = "white", high = "blue",name = "Number")+
  facet_grid(. ~ Status)+
  theme_classic()+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text=element_text(size=18,colour='black'),
        axis.title=element_text(size=20,colour='black'),
        strip.background=element_rect(colour="gray", fill="white"),
        strip.text = element_text(size=20,colour='black')
  )
p
plotfile='./fig/RNA_SNP.meancount.heatmap.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)
plotfile='./fig/RNA_SNP.meancount.heatmap.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)

# change  proportion heatmap ----------------------------------------------------
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
options(scipen = 9)
mydata<-read.delim(file ='merge_GATK_10sample_count_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5') ~ "SCNT-ABE"
)
mydata<- mydata %>%
  dplyr::group_by(base_change,Status) %>%
  dplyr::mutate(mean_proportion = round(mean(proportion),2))

# mean(mydata$proportion[mydata$base_change=='C>T' & mydata$Status=='SCNT-ABE'])
# mean(mydata$proportion[mydata$base_change=='C>T' & mydata$Status=='SCNT'])
# mean(mydata$proportion[mydata$base_change=='A>G' & mydata$Status=='SCNT-ABE'])
# mean(mydata$proportion[mydata$base_change=='A>G' & mydata$Status=='SCNT'])

mydata$From <- str_split_fixed(mydata$base_change,"[>]",n=2)[,1]
mydata$To <- str_split_fixed(mydata$base_change,"[>]",n=2)[,2]
mydata<-unique(mydata[,5:8])

colnames(mydata)
tmpdf<-data.frame("Status"=c(rep('SCNT',4),rep('SCNT-ABE',4)),
                  "mean_proportion"=c(rep(0,8)),
                  "From"=c('A','T','G','C','A','T','G','C'),
                  "To"=c('A','T','G','C','A','T','G','C'))
mydata<-rbind(as.data.frame(mydata),tmpdf)
mydata$Status<-factor(mydata$Status, levels=c('SCNT', 'SCNT-ABE'))

p<-ggplot(mydata, aes(x = To, y=From)) +
  geom_tile(aes(fill = mean_proportion),colour = "black") +
  geom_text(aes(label = mean_proportion),size=8,colour='black') +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_continuous(low = "white", high = "blue",name = "Proportion")+
  facet_grid(. ~ Status)+
  theme_classic()+
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text=element_text(size=18,colour='black'),
        axis.title=element_text(size=20,colour='black'),
        strip.background=element_rect(colour="gray", fill="white"),
        strip.text = element_text(size=20,colour='black')
  )
p
plotfile='./fig/RNA_SNP.mean_proportion.heatmap.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)
plotfile='./fig/RNA_SNP.mean_proportion.heatmap.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)


# #########boxplot --------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
options(scipen = 9)
mydata<-read.delim(file ='merge_GATK_10sample_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
mydata$number1<-mydata$number/1000

breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=base_change, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum+1000),breaks = breaks)+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
   stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum+500,hide.ns = TRUE)
# stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
#                    method = "wilcox.test",method.args =list(alternative = "greater"),size=5,label.y=maximum-1,hide.ns = TRUE)
p
plotfile='./fig/GATKbase_change_barplot_PC_NCofftargetv1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/GATKbase_change_barplot_PC_NCofftargetv1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)

pcnumber<-c(mydata$number[mydata$Status=='SCNT-ABE' & mydata$base_change =='A>G'])
ncnumber<-c(mydata$number[mydata$Status=='SCNT' & mydata$base_change =='A>G'])
pcmean<-mean(pcnumber)
ncmean<-mean(ncnumber)
pcse<-sd(pcnumber)/sqrt(length(pcnumber))
ncse<-sd(ncnumber)/sqrt(length(ncnumber))


######v2
breaks<-pretty(range(mydata$number1), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=base_change, y=number1, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs (X1000)")+
  scale_y_continuous(limits = c(0,maximum+2),breaks = breaks)+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum+1,hide.ns = TRUE)
# stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
#                    method = "wilcox.test",method.args =list(alternative = "greater"),size=5,label.y=maximum-1,hide.ns = TRUE)
p
plotfile='./fig/GATKbase_change_barplot_PC_NCofftargetv2.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/GATKbase_change_barplot_PC_NCofftargetv2.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)


###########add sum 
data1<- mydata %>%
  dplyr::group_by(Individual) %>%
  dplyr::summarise(number = sum(number))

data1$base_change=rep('ALL',nrow(data1))
data1$Status<-case_when(
  data1$Individual %in% c('NC1','NC2','NC3','NC4','NC5') ~ "SCNT",
  data1$Individual %in% c('PC1','PC2','PC3','PC4','PC5') ~ "SCNT-ABE"
)

mydata<-rbind(data1,mydata[,1:4])
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$base_change<-factor(mydata$base_change, levels=c('ALL',"A>G","C>T","A>C","A>T","C>G","C>A"))

############SD and mean
sem <- function(x) sqrt(var(x)/length(x))
tmpdata<-mydata[mydata$base_change=='A>G',]
(ABE_mean<-mean(tmpdata$number[tmpdata$Status=='SCNT-ABE'])) ##7888
(ABE_sem<-sem(tmpdata$number[tmpdata$Status=='SCNT-ABE'])) ##1672
(SCNT_meam<-mean(tmpdata$number[tmpdata$Status=='SCNT'])) ##2297
(SCNT_sem<-sem(tmpdata$number[tmpdata$Status=='SCNT']))   ##367.1

(ABE_sd<-sd(tmpdata$number[tmpdata$Status=='SCNT-ABE'])) ##3738
(SCNT_sd<-sd(tmpdata$number[tmpdata$Status=='SCNT']))   ##820.8
############end SD and mean

mydata$number1<-mydata$number/1000

breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=base_change, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum+1000),breaks = breaks)+
  scale_x_discrete(labels=c('All','A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum+500,hide.ns = TRUE)
# stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
#                    method = "wilcox.test",method.args =list(alternative = "greater"),size=5,label.y=maximum-1,hide.ns = TRUE)
p
plotfile='./fig/GATKbase_change_barplot_PC_NCofftarget_sumv1.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 6)
plotfile='./fig/GATKbase_change_barplot_PC_NCofftarget_sumv1.png'
ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 6)

######v2
breaks<-pretty(range(mydata$number1), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=base_change, y=number1, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs (X1000)")+
  scale_y_continuous(limits = c(0,maximum+2),breaks = breaks)+
  scale_x_discrete(labels=c('ALL','A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum+1,hide.ns = TRUE)
# stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
#                    method = "wilcox.test",method.args =list(alternative = "greater"),size=5,label.y=maximum-1,hide.ns = TRUE)
p
plotfile='./fig/GATKbase_change_barplot_PC_NCofftarget_sumv2.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 6)
plotfile='./fig/GATKbase_change_barplot_PC_NCofftarget_sumv2.png'
ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 6)

for (i in unique(mydata$base_change)) {
  x<-mydata$number1[mydata$base_change==i & mydata$Status=='SCNT-ABE']
  y<-mydata$number1[mydata$base_change==i & mydata$Status=='SCNT']
  # ttest<-t.test(x,y)
  wtest<-wilcox.test(x,y)
  # print(ttest)
  print(i)
  print(wtest$p.value)
}

#########################
#########boxplot proportion
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ="merge_GATK_10sample_count_merge_proportion.txt",header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))

p<-ggplot(mydata, aes(x=base_change, y=proportion*100, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  scale_y_continuous(limits = c(0,80),breaks = seq(0,80,10))+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=80, hide.ns = TRUE)
p
plotfile='./fig/GATKbarplot_PC_NC_proportion.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/GATKbarplot_PC_NC_proportion.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)


######PC
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_GATK_10sample_count_merge_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))

BEdata<-mydata[mydata$Status=='SCNT-ABE',]
BEdata$base_change<-factor(BEdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(BEdata, aes(x=base_change, y=proportion*100, color=base_change)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  scale_y_continuous(breaks = seq(0,100,20))+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/GATK_base_change_barplot_onlyPCofftarget.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)
plotfile='./fig/GATK_base_change_barplot_onlyPCofftarget.png'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)


######NC
CTdata<-mydata[mydata$Status=='SCNT',]
x<-CTdata$number[CTdata$base_change=="A>G" & CTdata$Status=='NC']
y<-CTdata$number[CTdata$base_change=="A>G" & CTdata$Status=='NC']
t.test(x,y)$p.value

CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(CTdata, aes(x=base_change, y=proportion*100, color=base_change)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  scale_y_continuous(limits = c(0,100),breaks = seq(0,100,20))+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/GATK_base_change_barplot_onlyNCofftarget.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)
plotfile='./fig/GATK_base_change_barplot_onlyNCofftarget.png'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)


#############number
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_GATK_10sample_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('NC1','NC2','NC3','NC4','NC5') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
mydata$number<-mydata$number/1000

BEdata<-mydata[mydata$Status=='SCNT-ABE',]
maximum<-max(BEdata$number)
maximum
breaks<-pretty(range(0,maximum*1.5), 8)
maximum<- breaks[length(breaks)]

BEdata$base_change<-factor(BEdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(BEdata, aes(x=base_change, y=number, color=base_change)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs (X1000)")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/GATK_count_barplot_onlyPCofftarget_number.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)
plotfile='./fig/GATK_count_barplot_onlyPCofftarget_number.png'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)

######NC
CTdata<-mydata[mydata$Status=='SCNT',]
CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))
maximum<-max(CTdata$number)
maximum
breaks<-pretty(range(0,maximum*2), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(CTdata, aes(x=base_change, y=number, color=base_change)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs (X1000)")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/GATK_count_barplot_onlyNCofftarget_number.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)
plotfile='./fig/GATK_count_barplot_onlyNCofftarget_number.png'
ggsave(plotfile, plot=p, dpi = 600,width = 6, height = 6)


#########jitterplot
rm(list=ls())
################CI
################How to Calculate Confidence Interval in R
#####https://datasharkie.com/how-to-calculate-confidence-interval-in-r/
###install.packages("Rmisc")
library(Rmisc)
mydata<-read.delim(file ='merge_GATK_10sample_Jitter.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ABE_data<-mydata[mydata$Status=='PC',]
CI(ABE_data$percent,    ci=0.95)
# upper  mean lower 
# 30.50 30.28 30.07 
summary(ABE_data$percent)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.0    17.6    21.9    30.3    31.2   100.0
min(ABE_data$percent) ##6

NC_data<-mydata[mydata$Status=='NC',]
CI(NC_data$percent,    ci=0.95)
# upper  mean lower 
# 42.82 42.24 41.66 
summary(NC_data$percent)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.67   18.60   25.71   42.24   62.50  100.00
min(NC_data$percent) ###6.67

all<-c(7784,13928,7684,3814,6230)
exon<-c(2245,5176,2777,1189,2228)
percnt<-sum(exon)/sum(all)
percnt1<-exon/all
############################V2#################
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_GATK_10sample_Jitter.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
Individual_label<-c("SCNT1","SCNT2","SCNT3","SCNT4","SCNT5","SCNT-ABE1","SCNT-ABE2","SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")
mydata<-mydata[order(mydata$Individual),]
breaks<-pretty(range(mydata$percent), 5)
maximum<- breaks[length(breaks)]
total_num<-table(mydata$Individual)
p<-ggplot(data = mydata, mapping = aes(x = Individual, y = percent,color=Status)) + 
  geom_jitter(width = 0.4,size = 0.001)+
  annotate(geom="text", x=c(1:10), y=103, label=total_num,
           color="black",size=6)+
  scale_y_continuous(breaks =breaks)+
  scale_x_discrete(labels=Individual_label)+
  labs(x = '',y = "RNA\nA-to-G editing (%)")+
  guides(color=FALSE)+
  scale_color_manual(values =c('#0571b0','#ca0020'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=15,color ="black"),axis.title=element_text(size=20,face="bold"),axis.text.x=element_text(angle=25,hjust=1,vjust=1))
p
plotfile='./fig/merge_GATK_10sample_Jitter.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
plotfile='./fig/merge_GATK_10sample_Jitter.png'
ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)

############chrs plot
chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y")
alldata<-read.delim(file ='merge_GATK_10sample_Jitter.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
for (i in c("PC1","PC2","PC3","PC4","PC5")) {
  mydata<-alldata[alldata$Individual==i & alldata$chr %in% chrs,]
  breaks<-pretty(range(mydata$percent), 5)
  maximum<- breaks[length(breaks)]
  mydata$chr<-factor(mydata$chr, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))
  p<-ggplot(data = mydata, mapping = aes(x = chr, y = percent,colour=chr)) + geom_jitter(width = 0.4,size = 0.5)+
    scale_y_continuous(limits = c(0,101),breaks =breaks)+
    labs(x = '',y = "RNA\nA-to-G editing (%)")+
    scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"))+
    theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'))+
    guides(color=FALSE)+
    theme(axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"))
  p
  plotfile=paste('./fig/GATK_',i,'_Jitter.pdf',sep = '')
  ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
  plotfile=paste('./fig/GATK_',i,'_Jitter.png',sep = '')
  ggsave(plotfile, plot=p, dpi = 600,width = 15, height = 10)
}


# # ggseqlogo######Load the required packages -----------------------------
# ggseqlogo######Load the required packages
# ggseqlogo######Load the required packages
rm(list=ls())
require(ggplot2)
require(ggseqlogo)
library(ggplot2)
library(ggpubr)

for (i in c('SCNT','ABE')) {
  fname=paste('all_',i,'_RNA_edit_siteseq.txt',sep='')
  mydata<-read.delim(file =fname,header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  colnames(mydata)<-c('pos','seq')
  p<-ggseqlogo(mydata$seq,method='p')
  plotfile=paste('./fig/all_',i,'_RNAedit_ggseqlogoAG.pdf',sep = '')
  ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
  plotfile=paste('./fig/all_',i,'_RNAedit_ggseqlogoAG.png',sep = '')
  ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
}

samplelist<-c("PC1","PC2","PC3","PC4","PC5","NC1","NC2","NC3","NC4","NC5")
for (i in samplelist) {
  fname=paste(i,'.RNA_edit_siteseq.txt',sep='')
  mydata<-read.delim(file =fname,header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
  colnames(mydata)<-c('pos','seq')
  p<-ggseqlogo(mydata$seq,method='p')
  plotfile=paste('./fig/',i,'.RNA_ggseqlogoAG.pdf',sep = '')
  ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
  plotfile=paste('./fig/',i,'.RNA_ggseqlogoAG.png',sep = '')
  ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
}

# ################gene expression -----------------------------------------
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(reshape2)
library(stringr)
options(scipen = 9)
mygene<-read.table(file ='A2Ggene.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mygene)<-c('gene','sample')
data_exp<-read.table(file ='/home/devdata/nyy/nyy_GFP/RNA_offtarget/count_TPM_ypp/TPM.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
rownames(data_exp)<-data_exp$Geneid
data_exp<-data_exp[,-1]
data_exp<-data_exp[rowSums(data_exp)>0,]
data_exp<-log2(data_exp+1)
data_exp<-na.omit(data_exp)
colnames(data_exp)<-str_split(colnames(data_exp),'[.]',simplify=T)[,2]

########PC1
set.seed(1234)
offtarget<-intersect(mygene$gene[mygene$sample=='PC1'],rownames(data_exp))
genes<-setdiff(rownames(data_exp),offtarget)
randoms<-sample(genes,length(offtarget))
PC1data<-data.frame(PC1=data_exp[offtarget,'PC1'], 
                    random_genes=data_exp[randoms,'PC1'], 
                    NC1=data_exp[offtarget,'NC1'], 
                    NC2=data_exp[offtarget,'NC2'], 
                    NC3=data_exp[offtarget,'NC3'], 
                    NC4=data_exp[offtarget,'NC4'], 
                    NC5=data_exp[offtarget,'NC5'],
                    rNC1=data_exp[randoms,'NC1'], 
                    rNC2=data_exp[randoms,'NC2'], 
                    rNC3=data_exp[randoms,'NC3'], 
                    rNC4=data_exp[randoms,'NC4'], 
                    rNC5=data_exp[randoms,'NC5'],stringsAsFactors = F)
PC1data$NC<-rowMeans(PC1data[,3:7])
PC1data$RNC<-rowMeans(PC1data[,8:12])
PC1data<-PC1data[,c('PC1','random_genes','NC','RNC')]

PC1data <- melt(PC1data,  variable.name = "sample", value.name = "TPM")

breaks<-pretty(range(PC1data$TPM,na.rm=T), 8)
maximum<- max(PC1data$TPM,na.rm = T)*1.5
my_comparisons <- list(c('PC1', 'random_genes'),c('PC1', "NC"),c('NC', "RNC"))
p<-ggplot(PC1data,aes(x=sample, y=TPM, fill=sample)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_jitter(size=0.4,position=position_jitter(0.1)) +
  # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "",y = "log2(TPM+1)")+
  scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('SCNT-ABE1\n(Offtarget genes)', 'SCNT-ABE1\n(Random genes)', 'All SCNT\n(Offtarget genes)', 'All SCNT\n(Random genes)'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  guides(fill=FALSE)+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p

plotfile='./fig/PC1_RNA_offtarget_expv2_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 9, height = 6)
plotfile='./fig/PC1_RNA_offtarget_expv2_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 9, height = 6)

########PC2
set.seed(1234)
offtarget<-intersect(mygene$gene[mygene$sample=='PC2'],rownames(data_exp))
genes<-setdiff(rownames(data_exp),offtarget)
randoms<-sample(genes,length(offtarget))
PC2data<-data.frame(PC2=data_exp[offtarget,'PC2'], 
                    random_genes=data_exp[randoms,'PC2'], 
                    NC1=data_exp[offtarget,'NC1'], 
                    NC2=data_exp[offtarget,'NC2'], 
                    NC3=data_exp[offtarget,'NC3'], 
                    NC4=data_exp[offtarget,'NC4'], 
                    NC5=data_exp[offtarget,'NC5'],
                    rNC1=data_exp[randoms,'NC1'], 
                    rNC2=data_exp[randoms,'NC2'], 
                    rNC3=data_exp[randoms,'NC3'], 
                    rNC4=data_exp[randoms,'NC4'], 
                    rNC5=data_exp[randoms,'NC5'],stringsAsFactors = F)
PC2data$NC<-rowMeans(PC2data[,3:7])
PC2data$RNC<-rowMeans(PC2data[,8:12])
PC2data<-PC2data[,c('PC2','random_genes','NC','RNC')]

PC2data <- melt(PC2data,  variable.name = "sample", value.name = "TPM")

breaks<-pretty(range(PC2data$TPM,na.rm=T), 8)
maximum<- max(PC2data$TPM,na.rm = T)*1.5
my_comparisons <- list(c('PC2', 'random_genes'),c('PC2', "NC"),c('NC', "RNC"))
p<-ggplot(PC2data,aes(x=sample, y=TPM, fill=sample)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_jitter(size=0.4,position=position_jitter(0.1)) +
  # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "",y = "log2(TPM+1)")+
  scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('SCNT-ABE1\n(Offtarget genes)', 'SCNT-ABE1\n(Random genes)', 'All SCNT\n(Offtarget genes)', 'All SCNT\n(Random genes)'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  guides(fill=FALSE)+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/PC2_RNA_offtarget_expv2_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC2_RNA_offtarget_expv2_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

########PC3
set.seed(1234)
offtarget<-intersect(mygene$gene[mygene$sample=='PC3'],rownames(data_exp))
genes<-setdiff(rownames(data_exp),offtarget)
randoms<-sample(genes,length(offtarget))
PC3data<-data.frame(PC3=data_exp[offtarget,'PC3'], 
                    random_genes=data_exp[randoms,'PC3'], 
                    NC1=data_exp[offtarget,'NC1'], 
                    NC2=data_exp[offtarget,'NC2'], 
                    NC3=data_exp[offtarget,'NC3'], 
                    NC4=data_exp[offtarget,'NC4'], 
                    NC5=data_exp[offtarget,'NC5'],
                    rNC1=data_exp[randoms,'NC1'], 
                    rNC2=data_exp[randoms,'NC2'], 
                    rNC3=data_exp[randoms,'NC3'], 
                    rNC4=data_exp[randoms,'NC4'], 
                    rNC5=data_exp[randoms,'NC5'],stringsAsFactors = F)
PC3data$NC<-rowMeans(PC3data[,3:7])
PC3data$RNC<-rowMeans(PC3data[,8:12])
PC3data<-PC3data[,c('PC3','random_genes','NC','RNC')]

PC3data <- melt(PC3data,  variable.name = "sample", value.name = "TPM")

breaks<-pretty(range(PC3data$TPM,na.rm=T), 8)
maximum<- max(PC3data$TPM,na.rm = T)*1.5
my_comparisons <- list(c('PC3', 'random_genes'),c('PC3', "NC"),c('NC', "RNC"))
p<-ggplot(PC3data,aes(x=sample, y=TPM, fill=sample)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_jitter(size=0.4,position=position_jitter(0.1)) +
  # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "",y = "log2(TPM+1)")+
  scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('SCNT-ABE1\n(Offtarget genes)', 'SCNT-ABE1\n(Random genes)', 'All SCNT\n(Offtarget genes)', 'All SCNT\n(Random genes)'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  guides(fill=FALSE)+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/PC3_RNA_offtarget_expv2_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC3_RNA_offtarget_expv2_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

########PC4
set.seed(1234)
offtarget<-intersect(mygene$gene[mygene$sample=='PC4'],rownames(data_exp))
genes<-setdiff(rownames(data_exp),offtarget)
randoms<-sample(genes,length(offtarget))
PC4data<-data.frame(PC4=data_exp[offtarget,'PC4'], 
                    random_genes=data_exp[randoms,'PC4'], 
                    NC1=data_exp[offtarget,'NC1'], 
                    NC2=data_exp[offtarget,'NC2'], 
                    NC3=data_exp[offtarget,'NC3'], 
                    NC4=data_exp[offtarget,'NC4'], 
                    NC5=data_exp[offtarget,'NC5'],
                    rNC1=data_exp[randoms,'NC1'], 
                    rNC2=data_exp[randoms,'NC2'], 
                    rNC3=data_exp[randoms,'NC3'], 
                    rNC4=data_exp[randoms,'NC4'], 
                    rNC5=data_exp[randoms,'NC5'],stringsAsFactors = F)
PC4data$NC<-rowMeans(PC4data[,3:7])
PC4data$RNC<-rowMeans(PC4data[,8:12])
PC4data<-PC4data[,c('PC4','random_genes','NC','RNC')]

PC4data <- melt(PC4data,  variable.name = "sample", value.name = "TPM")

breaks<-pretty(range(PC4data$TPM,na.rm=T), 8)
maximum<- max(PC4data$TPM,na.rm = T)*1.5
my_comparisons <- list(c('PC4', 'random_genes'),c('PC4', "NC"),c('NC', "RNC"))
p<-ggplot(PC4data,aes(x=sample, y=TPM, fill=sample)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_jitter(size=0.4,position=position_jitter(0.1)) +
  # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "",y = "log2(TPM+1)")+
  scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('SCNT-ABE1\n(Offtarget genes)', 'SCNT-ABE1\n(Random genes)', 'All SCNT\n(Offtarget genes)', 'All SCNT\n(Random genes)'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  guides(fill=FALSE)+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/PC4_RNA_offtarget_expv2_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC4_RNA_offtarget_expv2_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

########PC5
set.seed(1234)
offtarget<-intersect(mygene$gene[mygene$sample=='PC5'],rownames(data_exp))
genes<-setdiff(rownames(data_exp),offtarget)
randoms<-sample(genes,length(offtarget))
PC5data<-data.frame(PC5=data_exp[offtarget,'PC5'], 
                    random_genes=data_exp[randoms,'PC5'], 
                    NC1=data_exp[offtarget,'NC1'], 
                    NC2=data_exp[offtarget,'NC2'], 
                    NC3=data_exp[offtarget,'NC3'], 
                    NC4=data_exp[offtarget,'NC4'], 
                    NC5=data_exp[offtarget,'NC5'],
                    rNC1=data_exp[randoms,'NC1'], 
                    rNC2=data_exp[randoms,'NC2'], 
                    rNC3=data_exp[randoms,'NC3'], 
                    rNC4=data_exp[randoms,'NC4'], 
                    rNC5=data_exp[randoms,'NC5'],stringsAsFactors = F)
PC5data$NC<-rowMeans(PC5data[,3:7])
PC5data$RNC<-rowMeans(PC5data[,8:12])
PC5data<-PC5data[,c('PC5','random_genes','NC','RNC')]

PC5data <- melt(PC5data,  variable.name = "sample", value.name = "TPM")

breaks<-pretty(range(PC5data$TPM,na.rm=T), 8)
maximum<- max(PC5data$TPM,na.rm = T)*1.5
my_comparisons <- list(c('PC5', 'random_genes'),c('PC5', "NC"),c('NC', "RNC"))
p<-ggplot(PC5data,aes(x=sample, y=TPM, fill=sample)) +
  geom_boxplot(outlier.colour = NA) +
  # geom_jitter(size=0.4,position=position_jitter(0.1)) +
  # geom_dotplot(binaxis='y', stackdir='center',dotsize =0.4,position=position_dodge(0.75)) +
  labs(x = "",y = "log2(TPM+1)")+
  scale_y_continuous(limits = c(-0.5,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('SCNT-ABE1\n(Offtarget genes)', 'SCNT-ABE1\n(Random genes)', 'All SCNT\n(Offtarget genes)', 'All SCNT\n(Random genes)'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  guides(fill=FALSE)+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/PC5_RNA_offtarget_expv2_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/PC5_RNA_offtarget_expv2_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)




# ############gene and site overlap ------------------------------------------------
rm(list=ls())
library(ggplot2)
library(venn)
library(reshape2)
setwd("/home/devdata/nyy/nyy_GFP/RNA_offtarget/2021v3")
mydata<-read.table(file ='GATK_10sample_edit_site.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('sample','site')
# unique(mydata$sample)
A1 = unique(mydata$site[mydata$sample == 'PC1']) #把同 sample 的 gene 们归类到一个 list中，下略
B1 = unique(mydata$site[mydata$sample == 'PC2'])
C1 = unique(mydata$site[mydata$sample == 'PC3'])
D1 = unique(mydata$site[mydata$sample == 'PC4'])
E1 = unique(mydata$site[mydata$sample == 'PC5'])
intersect1<-Reduce(intersect,list(A1,B1,C1,D1,E1))
print (length(intersect1))
# print (c(length(A1),length(B1),length(C1),length(D1)))
A<-paste('SCNT-ABE1(',length(A1),')',sep='')
B<-paste('SCNT-ABE2(',length(B1),')',sep='')
C<-paste('SCNT-ABE3(',length(C1),')',sep='')
D<-paste('SCNT-ABE4(',length(D1),')',sep='')
E<-paste('SCNT-ABE5(',length(E1),')',sep='')
sname=paste(A,B,C,D,E,sep=',')
x = list(A1,B1,C1,D1,E1)   #将5组数据也就是5个 factor 放入一个 list 中
pdf('./fig/RNA_gene_overlap_sites_venn.pdf',width = 10, height = 8)
venn(x, snames = sname,ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8,ilcs=1.2)
dev.off()

########gene overlap
rm(list=ls())
library(ggplot2)
library(venn)
library(reshape2)
setwd("/home/devdata/nyy/nyy_GFP/RNA_offtarget/2021v3")
mydata<-read.table(file ='A2Ggene.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('gene','sample')
# unique(mydata$sample)

A1 = unique(mydata$gene[mydata$sample == 'PC1']) #把同 sample 的 gene 们归类到一个 list中，下略
B1 = unique(mydata$gene[mydata$sample == 'PC2'])
C1 = unique(mydata$gene[mydata$sample == 'PC3'])
D1 = unique(mydata$gene[mydata$sample == 'PC4'])
E1 = unique(mydata$gene[mydata$sample == 'PC5'])
intersect1<-Reduce(intersect,list(A1,B1,C1,D1,E1))
print (length(intersect1))
A<-paste('SCNT-ABE1(',length(A1),')',sep='')
B<-paste('SCNT-ABE2(',length(B1),')',sep='')
C<-paste('SCNT-ABE3(',length(C1),')',sep='')
D<-paste('SCNT-ABE4(',length(D1),')',sep='')
E<-paste('SCNT-ABE5(',length(E1),')',sep='')
sname=paste(A,B,C,D,E,sep=',')
x = list(A1,B1,C1,D1,E1)   #将5组数据也就是5个 factor 放入一个 list 中
pdf('./fig/RNA_gene_overlap_venn.pdf',width = 10, height = 8)
venn(x, snames = sname,ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 0.8,ilcs=1.2)
dev.off()

write.table(intersect1, file = "./RNA_gene_overlap.tab", quote = FALSE,sep="\t",row.names = F)


# 1.all gene PCA -----------------------------------------------------------------
rm(list = ls())
library(stringr)
edata<-read.delim(file = "/home/devdata/nyy/nyy_GFP/RNA_offtarget/count_TPM_ypp/TPM.tab",header = TRUE, sep = "\t",stringsAsFactors=FALSE)# 30770
sum(edata$Geneid == 'NA')
colnames(edata)<-c("Geneid","SCNT1","SCNT2","SCNT3","SCNT4","SCNT5","SCNT-ABE1","SCNT-ABE2",
                   "SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")
edata.matrix = as.matrix(edata[,-1])
ids = edata.matrix > 0  
index_0 = rowSums(ids) > 0
edata_0 = edata[index_0,] #24469

index_1 = rowSums(ids[,1:5]) > 0  & rowSums(ids[,6:10]) > 0
edata_1 = edata[index_1,] #21205

write.table(edata_0, file = "./edata_rm0.tab", quote = FALSE,sep="\t",row.names = TRUE)
write.table(edata_1, file = "./edata_rm1.tab", quote = FALSE,sep="\t",row.names = TRUE)

# PCA ---------------------------------------------------------------------
library(ggord)
library(ggplot2)

dif_data <- read.table(file = 'edata_rm0.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data = dif_data[,-1]
sampledata<-data.frame(SampleID=colnames(dif_data),
                       Group=c(rep('SCNT',5),rep('SCNT-ABE',5)))

pca_group=factor(c(rep('SCNT',5),rep('SCNT-ABE',5)),levels = c('SCNT-ABE','SCNT'))

edata.pca <- prcomp(t(dif_data), scale. = TRUE)
p <- ggord(edata.pca, grp_in = pca_group, arrow=0, vec_ext =0,txt=NULL,cols=c('red','blue'),
           # poly = FALSE,polylntyp='solid', # twodash, solid, longdash, dotted, dotdash, dashed, blank
           ellipse_pro = 0.8,alpha_el = 0.2,alpha = 1) + theme_bw() + 
        theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(color = 'black',size = 14),
        axis.title.x = element_text(color = 'black',size = 14),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x = element_text(color = 'black',size = 12),
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(color = 'black', size = 0.5),
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 10)
  ) + 
  annotate(geom = 'segment', y = Inf,  x = -Inf, yend = Inf, xend = Inf, color = 'black', size = 1)+ 
  annotate(geom = 'segment', y = -Inf, x = Inf, yend = Inf, xend = Inf, color = 'black', size = 1)
p
ggsave('./fig/SCNTvsABE_PCA.pdf', plot=p, dpi = 600,width = 6, height = 6)

# 2.all gene heatmap -------------------------------------------------------------
rm(list=ls())
library(pheatmap)
dif_data = read.table(file = 'edata_rm0.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data = dif_data[,-1]
colnames(dif_data)<-c("SCNT1","SCNT2","SCNT3","SCNT4","SCNT5","SCNT-ABE1","SCNT-ABE2",
                     "SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")
Z_score <- (dif_data - apply(dif_data, 1, mean)) / apply(dif_data, 1, sd)

annotation_col <- data.frame(Group = c(rep('SCNT',5),rep('ABE',5)))
rownames(annotation_col) <- colnames(dif_data)
ancols = c('blue','red')
names(ancols) = c("SCNT","ABE")
ann_colors <- list(Group=ancols)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = F,clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",
         show_rownames=F,annotation_col=annotation_col,annotation_colors=ann_colors,fontsize=9,
        filename = "./fig/gene_heatmap.pdf")

rm(list=ls())
library(pheatmap)
dif_data = read.table(file = 'edata_rm0.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data = dif_data[,-1]
colnames(dif_data)<-c("SCNT1","SCNT2","SCNT3","SCNT4","SCNT5","SCNT-ABE1","SCNT-ABE2",
                      "SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")
Z_score <- (dif_data - apply(dif_data, 1, mean)) / apply(dif_data, 1, sd)
####cluster gene
mat <- Z_score
d = dist(mat, method = 'euclidean')
tree = hclust(d, method = 'complete')
v = cutree(tree, 1)[tree$order]
gaps = which((v[-1] - v[-length(v)]) != 0)
gene.cluster <- as.data.frame(v)

dt <- Z_score[rownames(gene.cluster),]
annotation_col <- data.frame(row.names = colnames(dt),
                             Group = c(rep('SCNT',5),rep('ABE',5)))
ancols.Group <- c('#e41a1c','#377eb8')
names(ancols.Group) <- c('ABE','SCNT')
ann_colors <- list(Group = ancols.Group)
pheatmap(dt,border_color=NA,color = colorRampPalette(c("blue", "white", "red"))(50),
         scale = "none",
         cluster_rows = F,clustering_distance_rows = "euclidean",
         cluster_cols = T,clustering_distance_cols = "euclidean",
         show_rownames=F,annotation_col=annotation_col,annotation_colors=ann_colors,
         fontsize=7,width = 4, height = 4,
         filename = "./fig/all.gene_heatmap.pdf")

# 3. DEG analysis (limma) -------------------------------------------------
rm(list=ls())
library(limma)
dif_data = read.table(file = 'edata_rm0.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
rownames(dif_data) = dif_data[,'Geneid']
group_list=c(rep('SCNT',5),rep('SCNT-ABE',5))
group_list <- factor(group_list,levels = c("SCNT","SCNT-ABE"))

dat <- dif_data[,-1]
design=model.matrix(~factor(group_list))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)

deg=topTable(fit,coef=2,adjust='BH',number = Inf)
deg$gene<-rownames(deg)
head(deg) 
deg$change = as.factor(ifelse(deg$P.Value < 0.05,
                              ifelse( deg$logFC > 0 ,'UP','DOWN' ),
                              'NOT'))

deg$significance = as.factor(ifelse( deg$P.Value  <= 0.05,
                                     ifelse( deg$P.Value> 0.01, 'TRUE', 'M_TRUE' ),
                                     'FALSE' ))
table(deg$change)
# DOWN   NOT    UP 
# 320 23878   271
sum(deg$P.Value < 0.05)
# 591
table(deg$significance)
# FALSE M_TRUE   TRUE 
# 23878    103    488 
dif_data$gene<-rownames(dif_data)
deg_limma<-merge(dif_data,deg,by='gene')

ID_name<-read.table(file = "Macaca_mulatta_revised_ID_name.txt",header = TRUE, sep = "\t",stringsAsFactors=FALSE)# 30770
colnames(ID_name)<-c("Geneid","gene")
ID_name<-ID_name[ID_name$Geneid %in% deg_limma$Geneid,]
deg_limma<-merge(deg_limma,ID_name,by='Geneid',all = T)
colnames(deg_limma)[21]<-'gene_name'
write.table(deg_limma[,-2], file = "SCNT_deg_limma.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab

# exonic_variant_type plot ------------------------------------------------------
rm(list=ls())
library(dplyr)
setwd("/home/devdata/nyy/nyy_GFP/RNA_offtarget/2021v3")
fname=paste('PC1','_exonic_variant_infowithtype.txt',sep='')
df=read.delim(fname,header = F, stringsAsFactors = F)
colnames(df)<-c('SNV_type','gene','exon','ntalt','aaalt','rate')
df$sample<-'PC1'
all_ABE_df<-df

for (i in c('PC2','PC3','PC4','PC5')) {
  fname=paste(i,'_exonic_variant_infowithtype.txt',sep='')
  df=read.delim(fname,header = F, stringsAsFactors = F)
  colnames(df)<-c('SNV_type','gene','exon','ntalt','aaalt','rate')
  df$sample<-i
  all_ABE_df<-rbind(all_ABE_df,df)
}

mydata<-all_ABE_df[order(all_ABE_df$sample),]
mydata <- mydata %>% 
  mutate(SNV_type1 = case_when(
    .$SNV_type %in% c('stopgain')  ~ "stopgain",
    .$SNV_type %in% c('stoploss')  ~ "stoploss",
    .$SNV_type %in% c('nonsynonymous SNV')  ~ "nonsynonymous",
    .$SNV_type %in% c('synonymous SNV')  ~ "synonymous",
  )
)
mydata$SNV_type1<-factor(mydata$SNV_type1, levels=c('stopgain', 'stoploss','nonsynonymous','synonymous'))

prop_df<- mydata %>%
  dplyr::group_by(sample,SNV_type1) %>%
  dplyr::summarise(count = n() ) %>%
  dplyr::mutate( prop = count / sum(count) )

sample_label<-c("SCNT-ABE1","SCNT-ABE2","SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")
breaks<-pretty(range(mydata$rate), 5)
maximum<- breaks[length(breaks)]
total_num<-table(mydata$sample)
total_num<-paste(total_num,' sites',sep='')
p<-ggplot(data = mydata, mapping = aes(x = sample, y = rate,color=SNV_type1)) + 
  geom_jitter(width = 0.4,size = 1)+
  annotate(geom="text", x=c(1:5), y=103, label=total_num,
           color="black",size=6)+
  scale_y_continuous(breaks =breaks)+
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3'), name="")+
  scale_x_discrete(labels=sample_label)+
  labs(x = '',y = "RNA\nA-to-G editing (%)")+
  # scale_color_manual(values =c('#0571b0','#ca0020'))+
  theme_bw()+theme(legend.text=element_text(size=16),legend.title = element_blank(),legend.position = 'top')+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=15,color ="black"),axis.title=element_text(size=20,face="bold"),axis.text.x=element_text(angle=25,hjust=1,vjust=1))
p
ggsave(plot = p, './fig/exonic_varianttype_Editing_rate_Jitterplot.pdf', width = 8, height = 8,dpi = 600)

# fill histon plot ------------------------------------------------------
dt.freq <- as.data.frame(table(paste(mydata$sample,mydata$SNV_type1,sep = '_')))
dt.freq$sample <- str_split(dt.freq$Var1,'_',simplify = T)[,1]
dt.freq$type <- str_split(dt.freq$Var1,'_',simplify = T)[,2]
dt.freq$type <- factor(dt.freq$type,levels=c('stopgain', 'stoploss','nonsynonymous','synonymous'))

p <- ggplot(dt.freq, aes(x=sample, y=Freq, fill=type))+geom_col(position = 'fill') + 
  theme_bw() + labs(x = '', y = 'Percentage (%)') +
  scale_fill_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3'), name="") +
  scale_x_discrete(labels = sample_label) +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 12,vjust = -1.2),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x=element_text(angle=25,hjust=1,vjust=1,size = 10,color = 'black'),
        panel.grid = element_blank(),
        # strip.background = NULL, # 分面标题背景
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),
        legend.text = element_text(face = 'bold',color = 'black',size = 10))
p
ggsave(plot = p, './fig/exonic_varianttype_pct.pdf', width = 6, height = 5,dpi = 600)

# exonic_variant Census_cancer_gene----------------------------------------------------------
rm(list=ls())
library(dplyr)
setwd("/home/devdata/nyy/nyy_GFP/RNA_offtarget/2021v3")
ID_name<-read.table(file = "Macaca_mulatta_revised_ID_name.txt",header = TRUE, sep = "\t",stringsAsFactors=FALSE)# 30770
colnames(ID_name)<-c("gene","Gene.Symbol")
cancer_gene<-read.delim('Census_cancer_gene.tsv',header = T, stringsAsFactors = F)
cancer_gene <- cancer_gene %>% 
  mutate(generole = case_when(
    .$Role.in.Cancer %in% c('oncogene','oncogene, fusion')  ~ "oncogene",
    .$Role.in.Cancer %in% c("TSG, fusion",'TSG')  ~ "TSG",
    .$Role.in.Cancer %in% c("oncogene, TSG",'oncogene, TSG, fusion')  ~ "both",
    TRUE ~ "other"
  )
  )

table(cancer_gene$generole)
# both oncogene    other      TSG 
# 72      243      165      243 
cancer_gene1<-cancer_gene[,c("Gene.Symbol","generole")]
cancer_gene1<-merge(ID_name,cancer_gene1,by="Gene.Symbol", all = T)

fname=paste('PC1','_exonic_variant_info.txt',sep='')
df=read.delim(fname,header = F, stringsAsFactors = F)
colnames(df)<-c('gene','exon','ntalt','aaalt','rate')
df=merge(df,cancer_gene1,by='gene',all = F)
df$sample<-'PC1'
all_ABE_df<-df

for (i in c('PC2','PC3','PC4','PC5')) {
  fname=paste(i,'_exonic_variant_info.txt',sep='')
  df=read.delim(fname,header = F, stringsAsFactors = F)
  colnames(df)<-c('gene','exon','ntalt','aaalt','rate')
  df=merge(df,cancer_gene1,by='gene',all = F)
  df$sample<-i
  all_ABE_df<-rbind(all_ABE_df,df)
}
all_ABE_df<-all_ABE_df[! is.na(all_ABE_df$generole),]
all_ABE_df<-all_ABE_df[all_ABE_df$generole %in% c("oncogene","TSG","both"),]
write.table(all_ABE_df, file = "ABE_edit_rate.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab

common_site<-Reduce(intersect,list(all_ABE_df$ntalt[all_ABE_df$sample=='PC1'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC2'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC3'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC4'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC5'])
)

all_ABE_df_filter<-all_ABE_df[all_ABE_df$ntalt %in% common_site, ]
write.table(all_ABE_df_filter, file = "ABE_edit_rate_common.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab

# Editing rate Jitter point plot --------------------------------- --------
rm(list = ls())
library(dplyr)
library(stringr)
library(ggplot2)
mydata<- read.table(file = 'ABE_edit_rate.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
sample_label<-c("SCNT-ABE1","SCNT-ABE2","SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")
mydata$generole <- factor(mydata$generole,levels = c("both","oncogene","TSG"))
breaks<-pretty(range(mydata$rate), 5)
maximum<- breaks[length(breaks)]
total_num<-table(mydata$sample)
mean(total_num)
total_num<-paste(total_num,' sites',sep='')

p<-ggplot(data = mydata, mapping = aes(x = sample, y = rate,color=generole)) + 
  geom_jitter(width = 0.4,size = 1)+
  annotate(geom="text", x=c(1:5), y=103, label=total_num,
           color="black",size=6)+
  scale_y_continuous(breaks =breaks)+
  scale_x_discrete(labels=sample_label)+
  labs(x = '',y = "RNA\nA-to-G editing (%)")+
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a'), name="")+
  theme_bw()+theme(legend.text=element_text(size=16),legend.title = element_blank())+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=15,color ="black"),axis.title=element_text(size=20,face="bold"),axis.text.x=element_text(angle=25,hjust=1,vjust=1))
p
ggsave(plot = p, './fig/cancergene_Editing_rate_Jitterplot.pdf', width = 10, height = 10,dpi = 600)

#####fill histon plot
library(ggplot2)
dt <- mydata
dt.freq <- as.data.frame(table(paste(dt$sample,dt$generole,sep = ':')))
dt.freq$sample <- str_split(dt.freq$Var1,':',simplify = T)[,1]
dt.freq$type <- str_split(dt.freq$Var1,':',simplify = T)[,2]
dt.freq$type <- factor(dt.freq$type,levels = c("both","oncogene","TSG"))

p <- ggplot(dt.freq, aes(x=sample, y=Freq, fill=type))+geom_col(position = 'fill') + 
  theme_bw() + labs(x = '', y = 'Percentage (%)') +
  scale_fill_manual(values = c('#e41a1c','#377eb8','#4daf4a'), name="") +
  scale_x_discrete(labels = sample_label) +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 12,vjust = -1.2),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x=element_text(color = 'black',size = 10,angle=25,hjust=1,vjust=1),
        panel.grid = element_blank(),
        # strip.background = NULL, # 分面标题背景
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),
        legend.text = element_text(face = 'bold',color = 'black',size = 10))
p
ggsave(plot = p, './fig/ABE_cancergene_generole_pct.pdf', width = 5, height = 5,dpi = 600)

#####
rm(list = ls())
dt <- read.table(file = 'ABE_edit_rate_common.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
colnames(dt)
# [1] "gene"        "exon"        "ntalt"       "aaalt"       "rate"        "Gene.Symbol" "generole"    "sample"

#######CI,range,DEG check
library(Rmisc)
CI(dt$rate,ci=0.95)
# upper  mean lower 
# 66.22 57.58 48.93  
summary(dt$rate)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 16.1    35.7    54.1    57.6    83.8   100.0 
range(dt$rate)
# 16.15 100.00
DEG_df <- read.delim(file = 'SCNT_deg_limma.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
DRGlist<-DEG_df$gene_name[DEG_df$change !='NOT']
intersect(dt$Gene.Symbol,DRGlist) ###BMPR1A
#######end CI,range,DEG check

dt$generole <- factor(dt$generole,levels = c("both","oncogene","TSG"))
dt$ntalt <- gsub("[|]", "_", dt$ntalt)
dt$ntalt <- gsub(">", "/", dt$ntalt)
dt$ntalt <- paste('(chr',dt$ntalt,'_',sep = '')
dt$aaalt_new <- paste(substr(dt$aaalt,3,3),'/',
                      substr(dt$aaalt,nchar(dt$aaalt),nchar(dt$aaalt)),
                      ')',sep = '')

dt$name <- paste(dt$gene,dt$ntalt,dt$aaalt_new,sep = '')
dt$name <- factor(dt$name,levels = rev(unique(dt$name)))

dt$sample<-case_when(
  dt$sample =='PC1' ~ "SCNT-ABE1",
  dt$sample =='PC2' ~ "SCNT-ABE2",
  dt$sample =='PC3' ~ "SCNT-ABE3",
  dt$sample =='PC4' ~ "SCNT-ABE4",
  dt$sample =='PC5' ~ "SCNT-ABE5",
)

p <- ggplot(dt,aes(x=name,y=rate))+ coord_flip() + # 坐标轴翻转
  geom_bar(stat="identity",aes(fill=generole),width=0.6) + theme_bw() +
  labs(x = '', y = 'Editing rate (%)',title = '') +
  scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(0,100,20)), limits = c(0,110),expand = c(0,0)) +
  scale_fill_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628'), name="") +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 12),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 10),
        panel.grid = element_blank(),
        strip.background = NULL, # 分面标题背景
        strip.text = element_text(face = 'bold',color = 'black',size = 15),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 10))+
  facet_grid(. ~ sample,scales = 'free')
p
ggsave(plot = p, './fig/ABE_cancer_edit_rate_common_barplot_facet.pdf', width = 12, height = 4,dpi = 600)

####gene_logfc.plot
rm(list=ls())
library(dplyr)
TPM <- read.table(file = 'edata_rm0.tab', header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
rownames(TPM) <- TPM$Geneid
TPM <- TPM[,-1]
gene <- read.table(file = 'ABE_edit_rate_common.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
gene$generole <- factor(gene$generole,levels = c("both","oncogene","TSG"))
gene <- arrange(gene,generole)
gene <- gene[!duplicated(gene$Gene.Symbol),]
rownames(gene) <- gene$gene
dt <- TPM[rownames(gene),]
rownames(dt) <- gene$Gene.Symbol
annotation_col <- data.frame(row.names = colnames(dt),
                             Group = c(rep('SCNT',5),rep('ABE',5)))
ancols.Group <- c('#e41a1c','#377eb8')
names(ancols.Group) <- c('ABE','SCNT')
ann_colors <- list(Group = ancols.Group)
pheatmap(dt,border_color=NA,color = colorRampPalette(c("blue", "white", "red"))(50),
         scale = "row",
         cluster_rows = F,clustering_distance_rows = "euclidean",
         cluster_cols = T,clustering_distance_cols = "euclidean",
         show_rownames=T,annotation_col=annotation_col,annotation_colors=ann_colors,
         filename = "./fig/ABE_common_cancer_gene_heatmap.pdf",
         fontsize=7,width = 4.5, height = 5,)

dt$log2FC <- ''
for (i in 1:nrow(dt)) {
  dt$log2FC[i] <- log2(mean(as.numeric(dt[i,6:10]))/mean(as.numeric(dt[i,1:5])))
}
dt$gene <- rownames(dt)
dt$gene <- factor(dt$gene,levels = rev(dt$gene))
dt$log2FC <- as.numeric(dt$log2FC)
dt$change <- ifelse(dt$log2FC > 0 ,'UP','DOWN')
dt$change <- factor(dt$change,levels = c('UP','DOWN'))

library(ggplot2)
p <- ggplot(dt,aes(x=gene,y=log2FC,fill=change))+ coord_flip() + # 坐标轴翻转
  geom_bar(stat="identity",width=0.6) + theme_bw() +
  labs(x = '', y = 'log2FC',title = '') +
  # scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(0,100,20)), limits = c(0,110),expand = c(0,0)) +
  scale_fill_manual(values = c('#e41a1c','#377eb8'), name="") +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 6),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 6),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 6),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 6),
        # panel.grid = element_blank(),
        # strip.background = NULL, # 分面标题背景
        # strip.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 6))+
  geom_hline(yintercept = 0,lty=7,lwd=0.4,alpha=1)
p
ggsave(plot = p, './fig/ABE_common_cancer_gene_logfc.plot.pdf', width = 3.5, height = 4,dpi = 600)


# ####non coding RNA ------------------------------------------------------
rm(list=ls())
library(dplyr)
setwd("/home/devdata/nyy/nyy_GFP/RNA_offtarget/2021v3")
M_gene<-read.delim(file = 'Macaca_mulatta_gene_info.txt',header = F,sep = '\t',stringsAsFactors = F)
colnames(M_gene)<-c('Mgene','Mgenename','code','chr')
M_gene <- M_gene %>% 
  mutate(genecode = case_when(
    .$code %in% c("protein_coding")  ~ "protein_coding",
    .$code %in% c("lincRNA")  ~ "lincRNA",
    .$code %in% c("miRNA")  ~ "miRNA",
    TRUE ~ "other"
  ))

table(M_gene$genecode)
# lincRNA          miRNA          other protein_coding 
# 2938           2339           6010          21099 
all_ABE_df=read.delim('A2Ggene_freq.txt',header = F, stringsAsFactors = F)
colnames(all_ABE_df)<-c('sample','Mgene','rate','chr','pos','alt')
all_ABE_df=merge(all_ABE_df,M_gene,by='Mgene',all = F)

write.table(all_ABE_df, file = "ABE_noncode_edit_rate.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab
# table(all_ABE_df$genecode)
# lincRNA          miRNA protein_coding 
# 462              3          23934 

common_site<-Reduce(intersect,list(all_ABE_df$alt[all_ABE_df$sample=='PC1'],
                                   all_ABE_df$alt[all_ABE_df$sample=='PC2'],
                                   all_ABE_df$alt[all_ABE_df$sample=='PC3'],
                                   all_ABE_df$alt[all_ABE_df$sample=='PC4'],
                                   all_ABE_df$alt[all_ABE_df$sample=='PC5'])
)

all_ABE_noncode_df_filter<-all_ABE_df[all_ABE_df$alt %in% common_site, ]
write.table(all_ABE_noncode_df_filter, file = "Noncode_ABE_edit_rate_common.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab

###plot
mydata<-all_ABE_df
mydata<-mydata[mydata$genecode %in%c('protein_coding','lincRNA'),]
mydata$genecode <- factor(mydata$genecode,levels = c("protein_coding","lincRNA"))
mydata$sample<-factor(mydata$sample, levels=c('PC1','PC2','PC3','PC4','PC5'))
sample_label<-c("SCNT-ABE1","SCNT-ABE2","SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")

breaks<-pretty(range(mydata$rate), 5)
maximum<- breaks[length(breaks)]
total_num<-table(mydata$sample)
total_num<-paste(total_num,' sites',sep='')
p<-ggplot(data = mydata, mapping = aes(x = sample, y = rate,color=genecode)) + 
  geom_jitter(width = 0.4,size = 0.5)+
  annotate(geom="text", x=c(1:5), y=103, label=total_num,
           color="black",size=6)+
  scale_color_manual(values = c('#e41a1c','#377eb8'), name="") +
  scale_y_continuous(breaks =breaks)+
  scale_x_discrete(labels=sample_label)+
  labs(x = '',y = "RNA\nA-to-G editing (%)")+
  # scale_color_manual(values =c('#0571b0','#ca0020'))+
  theme_bw()+theme(legend.text=element_text(size=16),legend.title = element_blank())+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=15,color ="black"),axis.title=element_text(size=20,face="bold"),axis.text.x=element_text(angle=25,hjust=1,vjust=1))
p
ggsave(plot = p, './fig/Noncode_ABE_edit_rate_Jitterplot.pdf', width = 10, height = 6,dpi = 600)

#####fill histon plot
library(ggplot2)
dt <- mydata
dt.freq <- as.data.frame(table(paste(dt$sample,dt$genecode,sep = ':')))
dt.freq$sample <- str_split(dt.freq$Var1,':',simplify = T)[,1]
dt.freq$type <- str_split(dt.freq$Var1,':',simplify = T)[,2]
dt.freq <- dt.freq[dt.freq$type != 'miRNA',]

dt.freq$type <- factor(dt.freq$type,levels = c("protein_coding","lincRNA"))

p <- ggplot(dt.freq, aes(x=sample, y=Freq, fill=type))+geom_col(position = 'fill') + 
  theme_bw() + labs(x = '', y = 'Percentage (%)') +
  scale_fill_manual(values = c('#e41a1c','#377eb8'), name="") +
  scale_x_discrete(labels = sample_label) +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 12,vjust = -1.2),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x=element_text(color = 'black',size = 10,angle=25,hjust=1,vjust=1),
        panel.grid = element_blank(),
        # strip.background = NULL, # 分面标题背景
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),
        legend.text = element_text(face = 'bold',color = 'black',size = 10))
p
ggsave(plot = p, './fig/ABE_noncode_generole_pct.pdf', width = 5, height = 5,dpi = 600)

###start filter
rm(list=ls())
ABE_noncode<- read.table(file = 'ABE_noncode_edit_rate.tab', header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
H_gene<-read.delim(file = 'Homo_sapiens_gene_info.txt',header = F,sep = '\t',stringsAsFactors = F)
colnames(H_gene)<-c('Hgene','Hgenename','code','chr')
H_gene <- H_gene %>% 
  mutate(genecode = case_when(
    .$code %in% c("protein_coding")  ~ "protein_coding",
    .$code %in% c("lincRNA")  ~ "lincRNA",
    .$code %in% c("miRNA")  ~ "miRNA",
    TRUE ~ "other"
  ))

table(H_gene$genecode)
# lincRNA          miRNA          other protein_coding 
# 7668           4198          28983          19826
lncRNAdisease<-read.delim(file = 'experimentallncRNA-disease information.tab',header = T,sep = '\t',stringsAsFactors = F)
lncRNAdisease<-lncRNAdisease[lncRNAdisease$Species=='Homo sapiens',]
H_genev1<-H_gene[H_gene$Hgenename %in% unique(lncRNAdisease$ncRNA.Symbol),]
HM_gene<-read.delim(file = 'HM_gene_mart_export.txt',header = T,sep = '\t',stringsAsFactors = F)
colnames(HM_gene)<-c("Mgene","Hgene")
HM_gene_lnc<-merge(HM_gene,H_genev1,by="Hgene")

all_ABE_df=merge(ABE_noncode,HM_gene_lnc,by='Mgene',all = F)
all_ABE_df<-all_ABE_df[! is.na(all_ABE_df$genecode.x),]
# table(all_ABE_df$genecode.x)
# protein_coding 
# 323

# ####development gene ------------------------------------------- --------
rm(list=ls())
library(dplyr)
setwd("/home/devdata/nyy/nyy_GFP/RNA_offtarget/2021v3")
M_gene<-read.delim(file = 'Macaca_mulatta_gene_info.txt',header = F,sep = '\t',stringsAsFactors = F)
colnames(M_gene)<-c('Mgene','official_symbol','code','chr')
M_gene$gene<-M_gene$official_symbol
for (i in 1:nrow(M_gene)) {
  if (is.na(M_gene[i,'gene'])){
    M_gene[i,'gene']=M_gene[i,'Mgene']
  }
}

Blast_gene<-read.delim(file = 'Blast_markers.txt',header =T,sep = '\t',stringsAsFactors = F)
MBlast_gene<-merge(M_gene,Blast_gene,by='official_symbol')
MBlast_gene<-MBlast_gene[,c("official_symbol","Mgene","gene","cluster")]
dev_gene<-read.delim(file = 'ICSI.rmAbnormal.markers.txt',header =T,sep = '\t',stringsAsFactors = F)
dev_gene<-dev_gene[dev_gene$p_val_adj<0.01,]
Mdev_gene<-merge(M_gene,dev_gene,by='gene')
Mdev_gene<-Mdev_gene[,c("official_symbol","Mgene","gene","cluster")]
Mdev_gene<-rbind(Mdev_gene,MBlast_gene)
# table(Mdev_gene$cluster)
# 16C    1C    2C    4C    8C Blast   EPI     M   PrE    TE 
# 592   540  1693  1960   880   753   286   313   184   239 
gene_1C<-Mdev_gene$gene[Mdev_gene$cluster=='1C']
gene_2C<-Mdev_gene$gene[Mdev_gene$cluster=='2C']
gene_4C<-Mdev_gene$gene[Mdev_gene$cluster=='4C']
gene_8C<-Mdev_gene$gene[Mdev_gene$cluster=='8C']
gene_16C<-Mdev_gene$gene[Mdev_gene$cluster=='16C']
gene_morula<-Mdev_gene$gene[Mdev_gene$cluster=='M']
gene_TE<-Mdev_gene$gene[Mdev_gene$cluster=='TE']
gene_PrE<-Mdev_gene$gene[Mdev_gene$cluster=='PrE']
gene_EPI<-Mdev_gene$gene[Mdev_gene$cluster=='EPI']

gene_1Cdif<-Reduce(setdiff,list(gene_1C,gene_2C,gene_4C,gene_8C,gene_16C,
                                gene_morula,gene_TE,gene_PrE,gene_EPI)
)
gene_2Cdif<-Reduce(setdiff,list(gene_2C,gene_1C,gene_4C,gene_8C,gene_16C,
                                gene_morula,gene_TE,gene_PrE,gene_EPI)
)

gene_4Cdif<-Reduce(setdiff,list(gene_4C,gene_1C,gene_2C,gene_8C,gene_16C,
                                gene_morula,gene_TE,gene_PrE,gene_EPI)
)
gene_8Cdif<-Reduce(setdiff,list(gene_8C,gene_1C,gene_2C,gene_4C,gene_16C,
                                gene_morula,gene_TE,gene_PrE,gene_EPI)
)
gene_16Cdif<-Reduce(setdiff,list(gene_16C,gene_1C,gene_2C,gene_4C,gene_8C,
                                 gene_morula,gene_TE,gene_PrE,gene_EPI)
)
gene_moruladif<-Reduce(setdiff,list(gene_morula,gene_1C,gene_2C,gene_4C,
                                    gene_8C,gene_16C,gene_TE,gene_PrE,gene_EPI)
)
gene_TEdif<-Reduce(setdiff,list(gene_TE,gene_1C,gene_2C,gene_4C,gene_8C,
                                gene_16C,gene_morula,gene_PrE,gene_EPI)
)
gene_PrEdif<-Reduce(setdiff,list(gene_PrE,gene_1C,gene_2C,gene_4C,gene_8C,
                                gene_16C,gene_morula,gene_TE,gene_EPI)
)
gene_EPIdif<-Reduce(setdiff,list(gene_EPI,gene_1C,gene_2C,gene_4C,gene_8C,
                                gene_16C,gene_morula,gene_TE,gene_PrE)
)

Mdev_gene <- Mdev_gene %>%
  mutate(generole = case_when(
    .$gene %in% gene_1Cdif  ~ "1cell",
    .$gene %in% gene_2Cdif  ~ "2cell",
    .$gene %in% gene_4Cdif  ~ "4cell",
    .$gene %in% gene_8Cdif  ~ "8cell",
    .$gene %in% gene_16Cdif  ~ "16cell",
    .$gene %in% gene_moruladif  ~ "morula",
    .$gene %in% gene_TEdif  ~ "blastocyst-TE",
    .$gene %in% gene_PrEdif  ~ "blastocyst-PrE",
    .$gene %in% gene_EPIdif  ~ "blastocyst-EPI",
    TRUE ~ "other"
  ))
# table(Mdev_gene$generole)
# 16cell          1cell          2cell          4cell          8cell blastocyst-EPI blastocyst-PrE 
# 127             50            211            450            331            227            200 
# blastocyst-TE         morula          other 
# 246            162           5436

fname=paste('PC1','_exonic_variant_info.txt',sep='')
df=read.delim(fname,header = F, stringsAsFactors = F)
colnames(df)<-c('Mgene','exon','ntalt','aaalt','rate')
df=merge(df,Mdev_gene,by='Mgene',all = F)
df$sample<-'PC1'
all_ABE_df<-df

for (i in c('PC2','PC3','PC4','PC5')) {
  fname=paste(i,'_exonic_variant_info.txt',sep='')
  df=read.delim(fname,header = F, stringsAsFactors = F)
  colnames(df)<-c('Mgene','exon','ntalt','aaalt','rate')
  df=merge(df,Mdev_gene,by='Mgene',all = F)
  df$sample<-i
  all_ABE_df<-rbind(all_ABE_df,df)
}
all_ABE_df<-all_ABE_df[all_ABE_df$generole !='other',]
write.table(all_ABE_df, file = "ABE_dev_edit_rate.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab

common_site<-Reduce(intersect,list(all_ABE_df$ntalt[all_ABE_df$sample=='PC1'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC2'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC3'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC4'],
                                   all_ABE_df$ntalt[all_ABE_df$sample=='PC5'])
)

all_ABE_df_filter<-all_ABE_df[all_ABE_df$ntalt %in% common_site, ]
write.table(all_ABE_df_filter, file = "ABE_dev_edit_rate_common.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab

# Editing rate Jitter point plot --------------------------------- --------
rm(list = ls())
library(dplyr)
library(stringr)
library(ggplot2)
mydata<- read.table(file = 'ABE_dev_edit_rate.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
sample_label<-c("SCNT-ABE1","SCNT-ABE2","SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")
mydata$sample<-factor(mydata$sample, levels=paste('PC',1:5,sep = ''))
mydata$generole<- factor(mydata$generole,levels = c("1cell","2cell","4cell","8cell","16cell","morula","blastocyst-TE","blastocyst-EPI","blastocyst-PrE"))
breaks<-pretty(range(mydata$rate), 5)
maximum<- breaks[length(breaks)]
total_num<-table(mydata$sample)
total_num<-paste(total_num,' sites',sep='')
p<-ggplot(data = mydata, mapping = aes(x = sample, y = rate,color=generole)) + 
  geom_jitter(width = 0.4,size = 1)+
  annotate(geom="text", x=c(1:5), y=103, label=total_num,
           color="black",size=6)+
  scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'), name="") +
  scale_y_continuous(breaks =breaks)+
  scale_x_discrete(labels=sample_label)+
  labs(x = '',y = "RNA\nA-to-G editing (%)")+
  # scale_color_manual(values =c('#0571b0','#ca0020'))+
  theme_bw()+theme(legend.text=element_text(size=16),legend.title = element_blank())+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=15,color ="black"),axis.title=element_text(size=20,face="bold"),axis.text.x=element_text(angle=25,hjust=1,vjust=1))
p
ggsave(plot = p, './fig/devgene_Editing_rate_Jitterplot.pdf', width = 10, height = 6,dpi = 600)

# fill Stacked bars plot ------------------------------------------------------
rm(list=ls())
library(ggplot2)
dt <- read.table(file = 'ABE_dev_edit_rate.tab', header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dt.freq <- as.data.frame(table(paste(dt$sample,dt$generole,sep = '_')))
dt.freq$sample <- str_split(dt.freq$Var1,'_',simplify = T)[,1]
dt.freq$type <- str_split(dt.freq$Var1,'_',simplify = T)[,2]
dt.freq$type <- factor(dt.freq$type,levels = c("1cell","2cell","4cell","8cell","16cell","morula","blastocyst-TE","blastocyst-EPI","blastocyst-PrE"))
sample_label<-c("SCNT-ABE1","SCNT-ABE2","SCNT-ABE3","SCNT-ABE4","SCNT-ABE5")


########mean_proportion for each stage key genes
prop_df<- dt %>%
  dplyr::group_by(sample,generole) %>%
  dplyr::summarise(count = n() ) %>%
  dplyr::mutate( prop = count / sum(count) )

prop_df<- prop_df %>%
  dplyr::group_by(generole) %>%
  dplyr::mutate(mean_proportion = round(mean(prop),2))
unique(prop_df[,c(2,5)])
########end mean_proportion for each stage key genes

p <- ggplot(dt.freq, aes(x=sample, y=Freq, fill=type))+geom_col(position = 'fill') + 
  theme_bw() + labs(x = '', y = 'Percentage (%)') +
  scale_fill_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'), name="") +
  scale_x_discrete(labels = sample_label) +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 12,vjust = -1.2),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x=element_text(angle=20,hjust=1,vjust=1,color = 'black',size = 10),
        panel.grid = element_blank(),
        # strip.background = NULL, # 分面标题背景
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),
        legend.text = element_text(face = 'bold',color = 'black',size = 10))
p
ggsave(plot = p, './fig/ABE_dev_edit_rate.generole_pct.pdf', width = 5, height = 5,dpi = 600)

###3.ABE_dev_edit_rate_common barplot_facet -------------------------------------------------
rm(list = ls())
library(dplyr)
library(stringr)
library(ggplot2)
dt <- read.table(file = 'ABE_dev_edit_rate_common.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
colnames(dt)

#######CI,range,DEG check
library(Rmisc)
CI(dt$rate,ci=0.95)
# upper  mean lower 
# 61.25 57.47 53.68  
summary(dt$rate)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14.9    38.1    57.6    57.5    74.3   100.0
range(dt$rate)
DEG_df <- read.delim(file = 'SCNT_deg_limma.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
DRGlist<-DEG_df$gene_name[DEG_df$change !='NOT']
intersect(dt$gene,DRGlist) ##FADS1
intersect(dt$official_symbol,DRGlist) ###"FADS1"
#######end CI,range,DEG check

dt$ntalt <- gsub("[|]", "_", dt$ntalt)
dt$ntalt <- gsub(">", "/", dt$ntalt)
dt$ntalt <- paste('(chr',dt$ntalt,'_',sep = '')
dt$aaalt_new <- paste(substr(dt$aaalt,3,3),'/',
                      substr(dt$aaalt,nchar(dt$aaalt),nchar(dt$aaalt)),
                      ')',sep = '')

dt$name <- paste(dt$gene,dt$ntalt,dt$aaalt_new,sep = '')
colnames(dt)
dt$generole <- factor(dt$generole,levels = c("1cell","2cell","4cell","8cell","16cell","morula","blastocyst-TE","blastocyst-EPI","blastocyst-PrE"))
dt <- arrange(dt,generole)
dt$name <- factor(dt$name,levels = rev(unique(dt$name)))

dt$sample<-case_when(
  dt$sample =='PC1' ~ "SCNT-ABE1",
  dt$sample =='PC2' ~ "SCNT-ABE2",
  dt$sample =='PC3' ~ "SCNT-ABE3",
  dt$sample =='PC4' ~ "SCNT-ABE4",
  dt$sample =='PC5' ~ "SCNT-ABE5",
)

p <- ggplot(dt,aes(x=name,y=rate))+ coord_flip() + # 坐标轴翻转
  geom_bar(stat="identity",aes(fill=generole),width=0.6) + theme_bw() +
  labs(x = '', y = 'Editing rate (%)',title = '') +
  scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(0,100,20)), limits = c(0,110),expand = c(0,0)) +
  scale_fill_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'), name="") +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 12),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 12),
        panel.grid = element_blank(),
        strip.background = NULL, # 分面标题背景
        strip.text = element_text(face = 'bold',color = 'black',size = 15),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 10))+
  facet_grid(. ~ sample,scales = 'free')
p
ggsave(plot = p, './fig/ABE_dev_edit_rate_common_barplot_facet.pdf', width = 15, height = 8,dpi = 600)

# development common gene FC 
rm(list=ls())
library(dplyr)
TPM <- read.table(file = 'edata_rm0.tab', header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
rownames(TPM) <- TPM$Geneid
TPM <- TPM[,-1]
gene <- read.table(file = 'ABE_dev_edit_rate_common.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
gene$generole <- factor(gene$generole,levels =c("1cell","2cell","4cell","8cell","16cell","morula","blastocyst-TE","blastocyst-EPI","blastocyst-PrE"))
gene <- arrange(gene,generole)
gene <- gene[!duplicated(gene$Mgene),]
rownames(gene) <- gene$Mgene
dt <- TPM[rownames(gene),]
rownames(dt) <- gene$gene
annotation_col <- data.frame(row.names = colnames(dt),
                             Group = c(rep('SCNT',5),rep('ABE',5)))
ancols.Group <- c('#e41a1c','#377eb8')
names(ancols.Group) <- c('ABE','SCNT')
ann_colors <- list(Group = ancols.Group)
pheatmap(dt,border_color=NA,color = colorRampPalette(c("blue", "white", "red"))(50),
         scale = "row",
         cluster_rows = F,clustering_distance_rows = "euclidean",
         cluster_cols = T,clustering_distance_cols = "euclidean",
         show_rownames=T,annotation_col=annotation_col,annotation_colors=ann_colors,
         filename = "./fig/ABE_dev_gene_heatmap.pdf",
         fontsize=7,width = 4.5, height = 5,)

dt$log2FC <- ''
for (i in 1:nrow(dt)) {
  dt$log2FC[i] <- log2(mean(as.numeric(dt[i,6:10]))/mean(as.numeric(dt[i,1:5])))
}
dt$gene <- rownames(dt)
dt$gene <- factor(dt$gene,levels = rev(dt$gene))
dt$log2FC <- as.numeric(dt$log2FC)
dt$change <- ifelse(dt$log2FC > 0 ,'UP','DOWN')
dt$change <- factor(dt$change,levels = c('UP','DOWN'))

library(ggplot2)
p <- ggplot(dt,aes(x=gene,y=log2FC,fill=change))+ coord_flip() + # 坐标轴翻转
  geom_bar(stat="identity",width=0.6) + theme_bw() +
  labs(x = '', y = 'log2FC',title = '') +
  # scale_y_continuous(breaks = seq(0,100,20),labels = as.character(seq(0,100,20)), limits = c(0,110),expand = c(0,0)) +
  scale_fill_manual(values = c('#e41a1c','#377eb8'), name="") +
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 6),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 6),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 6),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 6),
        # panel.grid = element_blank(),
        # strip.background = NULL, # 分面标题背景
        # strip.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 6))+
  geom_hline(yintercept = 0,lty=7,lwd=0.4,alpha=1)
p
ggsave(plot = p, './fig/ABE_dev_gene_logfc.plot.pdf', width = 3.5, height = 4,dpi = 600)

# 5.A2G gene GO/KEGG --------------------------------------------------
rm(list=ls())
library(clusterProfiler)
gene.dt <- read.table(file = 'A2Ggene.txt',sep='\t',header=F, stringsAsFactors=FALSE)
table(gene.dt$V2)
# PC1  PC2  PC3  PC4  PC5 
# 2749 4290 2961 1843 2550

enrich.res <- data.frame(ID='ID',Description='Description',GeneRatio='GeneRatio',BgRatio='BgRatio',pvalue='pvalue',
                         p.adjust='p.adjust',qvalue='qvalue',geneID='geneID',Count='Count',Type='Type',Group='Group',
                         stringsAsFactors = F)
enrich.res <- enrich.res[-1,]

load(file = "/nas01/Genome/enrichment_annotation/Rhesus_enrichment.RData")
term2gene.kegg <- kegg[,c(1,4)]
term2name.kegg <- kegg[,c(1,2)]

term2gene.bp <- gobp[,c(1,4)]
term2name.bp <- gobp[,c(1,2)]

term2gene.mf <- gomf[,c(1,4)]
term2name.mf <- gomf[,c(1,2)]

term2gene.cc <- gocc[,c(1,4)]
term2name.cc <- gocc[,c(1,2)]

for (i in unique(gene.dt$V2)) {
  gene <- gene.dt$V1[gene.dt$V2 == i]
  tag <- i
  
  if (length(gene) > 0) {
    # KEGG
    ekegg <- enricher(
      gene = gene,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 1,
      maxGSSize = 500,
      qvalueCutoff = 1,
      TERM2GENE = term2gene.kegg,
      TERM2NAME = term2name.kegg
    )
    
    if (length(ekegg) != 0) {
      ekegg_result<-as.data.frame(ekegg@result)
      ekegg_result$Type <- 'KEGG'
      ekegg_result$Group <- tag
      enrich.res <- rbind(enrich.res,ekegg_result)
    }
    
    # GO BP
    ego_BP <- enricher(
      gene = gene,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 1,
      maxGSSize = 500,
      qvalueCutoff = 1,
      TERM2GENE = term2gene.bp,
      TERM2NAME = term2name.bp
    )
    if (length(ego_BP) != 0) {
      ego_BP_result<-as.data.frame(ego_BP@result)
      ego_BP_result$Type <- 'GO_BP'
      ego_BP_result$Group <- tag
      enrich.res <- rbind(enrich.res,ego_BP_result)
    }
    
    # GO MF
    ego_MF <- enricher(
      gene = gene,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 1,
      maxGSSize = 500,
      qvalueCutoff = 1,
      TERM2GENE = term2gene.mf,
      TERM2NAME = term2name.mf
    )
    if (length(ego_MF) != 0) {
      ego_MF_result<-as.data.frame(ego_MF@result)
      ego_MF_result$Type <- 'GO_MF'
      ego_MF_result$Group <- tag
      enrich.res <- rbind(enrich.res,ego_MF_result)
    }
    
    # GO CC
    ego_CC <- enricher(
      gene = gene,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 1,
      maxGSSize = 500,
      qvalueCutoff = 1,
      TERM2GENE = term2gene.cc,
      TERM2NAME = term2name.cc
    )
    if (length(ego_CC) != 0) {
      ego_CC_result<-as.data.frame(ego_CC@result)
      ego_CC_result$Type <- 'GO_CC'
      ego_CC_result$Group <- tag
      enrich.res <- rbind(enrich.res,ego_CC_result)
    }
  }
}

write.table(enrich.res, file = 'A2Ggene.gene.enrich.res.txt', quote = FALSE,sep="\t",row.names = FALSE)

# # A2Ggene GO BP heatmap  ------------------------------------------------
rm(list = ls())   
library(ggplot2)
library(dplyr)
library(reshape2)
dt.all <- read.delim(file = 'A2Ggene.gene.enrich.res.txt',sep='\t',header=T, stringsAsFactors=FALSE)
dt.sub <- dt.all[dt.all$Type == 'GO_BP' & dt.all$pvalue < 0.05,c(2,5,9,11)] %>% arrange(Group,pvalue)

table(dt.sub$Group)
# PC1 PC2 PC3 PC4 PC5 
# 188 189 216 179 183

# get top term
n.top <- 10
for (i in unique(dt.sub$Group)) {
  tmp <- dt.sub[dt.sub$Group == i,]
  if (nrow(tmp) >= n.top) {
    tmp.top <- tmp[1:n.top,]
  } else {
    tmp.top <- tmp
  }
  
  if (!exists('dt.top')) {
    dt.top <- tmp.top
  } else {
    dt.top <- rbind(dt.top,tmp.top)
  }
}

# 补齐数据
dt.top.all <- dt.sub[dt.sub$Description %in% dt.top$Description,]
for (i in unique(dt.top$Description)) {
  tmp <- dt.top.all[dt.top.all$Description == i,]
  if (nrow(tmp) != 5) {
    n.cluster <- setdiff(unique(dt.top$Group),unique(tmp$Group))
    
    tmp.add <-dt.top.all[1:length(n.cluster),]
    tmp.add$Group <- n.cluster
    tmp.add$Description <- unique(tmp$Description)
    tmp.add$pvalue <- NA
    tmp.add$Count <- NA
    
    tmp <- rbind(tmp,tmp.add)
  }
  
  if (!exists('dt.top.all.add')) {
    dt.top.all.add <- tmp
  } else {
    dt.top.all.add <- rbind(dt.top.all.add,tmp)
  }
}

dt.shape <- melt(dt.top.all.add)
dt.shape$Group <- factor(dt.shape$Group,levels = unique(dt.top$Group))
dt.shape$tag <- paste(dt.shape$Group,dt.shape$Description,sep = '_')

dt.top.all.order <- dt.top.all[dt.top.all$pvalue <= 0.05,]
dt.shape$Description <- factor(dt.shape$Description,levels = rev(unique(dt.top.all.order$Description)))

dt.shape.count <- dt.shape[dt.shape$variable == 'Count',]
dt.shape.p <- dt.shape[dt.shape$variable == 'pvalue',]

dt.shape.p$value2 <- ifelse(dt.shape.p$value > 0.05,NA,-log10(dt.shape.p$value))
dt.shape.count$value2 <- ifelse(dt.shape.count$tag %in% dt.shape.p$tag[is.na(dt.shape.p$value2)],
                                NA,dt.shape.count$value)

p <- ggplot(dt.shape.p, aes(x=Group, y= Description)) + 
  geom_tile(aes(fill = value2),colour = "white") + 
  scale_fill_gradient2(name="-log10(pvalue)", guide = guide_colorbar(reverse = F),
                       midpoint = max(dt.shape.p$value2[!is.na(dt.shape.p$value2)])/2,
                       low = "#fee5d9",mid = "#fb6a4a",high = "#99000d",
                       n.breaks = 7, na.value = "grey80") +
  coord_fixed(ratio=0.25) + theme_minimal() +
  scale_x_discrete(labels=c('SCNT-ABE1','SCNT-ABE2','SCNT-ABE3','SCNT-ABE4','SCNT-ABE5'))+
  theme(axis.text.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.text.x=element_text(angle=20,hjust=1,vjust=1,color = 'black',size = 12),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.title = element_text(size = 8))+ 
  labs(x = "", y = "", title = "") +
  geom_point(data = dt.shape.count, aes(x=Group, y= Description,size=value2,color='#005a32')) +
  scale_color_manual(name="Count",values = c('#41ab5d')) +
  scale_size_continuous(name="Count", range=c(0.5,5), #guide = guide_legend(reverse = T),
                        breaks = round(seq(min(dt.shape.count$value2[!is.na(dt.shape.count$value2)]),
                                           max(dt.shape.count$value2[!is.na(dt.shape.count$value2)]), length.out = 6),digits = 0))

p
ggsave(plot = p, './fig/A2Ggene.gene_GO_BP.plot.pdf', width = 10, height = 10,dpi = 600)

# KEGG heatmap
rm(list = ls())   
library(ggplot2)
library(dplyr)
library(reshape2)

dt.all <- read.delim(file = 'A2Ggene.gene.enrich.res.txt',sep='\t',header=T, stringsAsFactors=FALSE)
dt.sub <- dt.all[dt.all$Type == 'KEGG' & dt.all$pvalue < 0.05,c(2,5,9,11)] %>% arrange(Group,pvalue)

table(dt.sub$Group)
# PC1 PC2 PC3 PC4 PC5 
# 28  47  30  19  50 

# get top term
n.top <- 10
for (i in unique(dt.sub$Group)) {
  tmp <- dt.sub[dt.sub$Group == i,]
  if (nrow(tmp) >= n.top) {
    tmp.top <- tmp[1:n.top,]
  } else {
    tmp.top <- tmp
  }
  
  if (!exists('dt.top')) {
    dt.top <- tmp.top
  } else {
    dt.top <- rbind(dt.top,tmp.top)
  }
}

# 补齐数据
dt.top.all <- dt.sub[dt.sub$Description %in% dt.top$Description,]
for (i in unique(dt.top$Description)) {
  tmp <- dt.top.all[dt.top.all$Description == i,]
  if (nrow(tmp) != 5) {
    n.cluster <- setdiff(unique(dt.top$Group),unique(tmp$Group))
    
    tmp.add <-dt.top.all[1:length(n.cluster),]
    tmp.add$Group <- n.cluster
    tmp.add$Description <- unique(tmp$Description)
    tmp.add$pvalue <- NA
    tmp.add$Count <- NA
    
    tmp <- rbind(tmp,tmp.add)
  }
  
  if (!exists('dt.top.all.add')) {
    dt.top.all.add <- tmp
  } else {
    dt.top.all.add <- rbind(dt.top.all.add,tmp)
  }
}

dt.shape <- melt(dt.top.all.add)
dt.shape$Group <- factor(dt.shape$Group,levels = unique(dt.top$Group))
dt.shape$tag <- paste(dt.shape$Group,dt.shape$Description,sep = '_')

dt.top.all.order <- dt.top.all[dt.top.all$pvalue <= 0.05,]
dt.shape$Description <- factor(dt.shape$Description,levels = rev(unique(dt.top.all.order$Description)))

dt.shape.count <- dt.shape[dt.shape$variable == 'Count',]
dt.shape.p <- dt.shape[dt.shape$variable == 'pvalue',]

dt.shape.p$value2 <- ifelse(dt.shape.p$value > 0.05,NA,-log10(dt.shape.p$value))
dt.shape.count$value2 <- ifelse(dt.shape.count$tag %in% dt.shape.p$tag[is.na(dt.shape.p$value2)],
                                NA,dt.shape.count$value)

p <- ggplot(dt.shape.p, aes(x=Group, y= Description)) + 
  geom_tile(aes(fill = value2),colour = "white") + 
  scale_fill_gradient2(name="-log10(pvalue)", guide = guide_colorbar(reverse = F),
                       midpoint = max(dt.shape.p$value2[!is.na(dt.shape.p$value2)])/2,
                       low = "#fee5d9",mid = "#fb6a4a",high = "#99000d",
                       n.breaks = 5, na.value = "grey80") +
  coord_fixed(ratio=0.25) + theme_minimal() +
  scale_x_discrete(labels=c('SCNT-ABE1','SCNT-ABE2','SCNT-ABE3','SCNT-ABE4','SCNT-ABE5'))+
  theme(axis.text.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.text.x=element_text(angle=20,hjust=1,vjust=1,color = 'black',size = 12),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.title = element_text(size = 8))+ 
  labs(x = "", y = "", title = "") +
  geom_point(data = dt.shape.count, aes(x=Group, y= Description,size=value2,color='#005a32')) +
  scale_color_manual(name="Count",values = c('#41ab5d')) +
  scale_size_continuous(name="Count", range=c(0.5,5), #guide = guide_legend(reverse = T),
                        breaks = round(seq(min(dt.shape.count$value2[!is.na(dt.shape.count$value2)]),
                                           max(dt.shape.count$value2[!is.na(dt.shape.count$value2)]), length.out = 6),digits = 0))
p
ggsave(plot = p, './fig/A2Ggene.gene_KEGG.plot.pdf', width = 10, height = 10,dpi = 600)

