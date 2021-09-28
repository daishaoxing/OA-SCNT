# raw count compare different condition ---------------------------------------------
rm(list=ls())
metadata<-read.delim(file ='sample_info.txt',header =T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-read.delim(file ='GATK_DP40_SNP_count.txt',header =T, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('Individual','Number1','Number2')
mydata<-merge(mydata,metadata,by='Individual',all = T)

allsample<-c("08431GFP","08431neg","08431pos","08431WT","090202body","090202head","cloning-JR","cloning-PF","cloning-W",
             "NC1","NC2","NC3","NC4","NC5","NC6","NC7","NC8","NC9","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")

###Donor cell vs  SCNT 
##The samples of donor cells (Promega #A1125) vs. SCNT embryos (QIAGEN, #150343)
library(ggplot2)
library(ggpubr)
options(scipen = 9)
samplelist<-c("08431GFP","08431neg","08431pos","08431WT","NC1","NC2","NC3","NC4","NC5","NC6","NC7","NC8","NC9")
mydata1<-mydata[mydata$Individual %in% samplelist,]
# mydata1<-mydata
mydata1$Number1<-mydata1$Number1/1000000
breaks<-pretty(range(mydata1$Number1), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=Status1, y=Number1, color=Status1)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=2)+
  labs(x = "",y = "Number of SNVs (x10e6)")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(aes(group = PCR,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/GATK_raw_count_dp40_barplot_SCNTvsdonor_cells_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/GATK_raw_count_dp40_barplot_SCNTvsdonor_cells_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)


##The samples of embryos (Promega #A1125)  vs. tissues. (QIAGEN, #150343).
library(ggplot2)
library(ggpubr)
options(scipen = 9)
samplelist<-c("cloning-JR","cloning-PF","cloning-W","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8",
              "090202body","090202head","NC1","NC2","NC3","NC4","NC5","NC6","NC7","NC8","NC9")
mydata1<-mydata[mydata$Individual %in% samplelist,]
mydata1$Number1<-mydata1$Number1/1000000
breaks<-pretty(range(mydata1$Number1), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=Dev, y=Number1, color=Dev)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=2)+
  labs(x = "",y = "Number of SNVs (x10e6)")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(aes(group = Dev,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/GATK_raw_count_dp40_barplot_Development_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/GATK_raw_count_dp40_barplot_Development_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

####Status
library(ggplot2)
library(ggpubr)
options(scipen = 9)
samplelist<-c("NC1","NC2","NC3","NC4","NC5","NC6","NC7","NC8",
              "NC9","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8")
mydata1<-mydata[mydata$Individual %in% samplelist,]
mydata1$Number1<-mydata1$Number1/1000000
breaks<-pretty(range(mydata1$Number1), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata1, aes(x=Status, y=Number1, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs (x10e6)")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/GATK_raw_count_dp40_barplot_Status_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/GATK_raw_count_dp40_barplot_Status_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)


# on-target editing from the WGS results ----------------------------------
rm(list=ls())
library(ggplot2)
library(reshape2)
samplelist<-c('SCNT1','SCNT2','SCNT3','SCNT4','SCNT5','SCNT6','SCNT7','SCNT8','SCNT9',
              'SCNT10-Body','SCNT11-Head',
              'SCNT-ABE1','SCNT-ABE2','SCNT-ABE3','SCNT-ABE4','SCNT-ABE5','SCNT-ABE6','SCNT-ABE7','SCNT-ABE8',
              'SCNT-ABE9-Muscle','SCNT-ABE10-Skin','SCNT-ABE11-Stomach')
editrate<-c(rep(0,11),72,90,100,94,100,95,100,100,100,100,100)
mydata<-data.frame(sample=samplelist,editrate=editrate)

breaks<-pretty(range(-1,mydata$editrate), 6)
maximum<- breaks[length(breaks)]
mydata$sample<-factor(mydata$sample, levels=samplelist)
rates<-mydata$editrate
# names(rates)<-mydata$sample
rates<-paste(rates,'%',sep='')
p<-ggplot(data = mydata, mapping = aes(x = sample, y = editrate,fill = sample)) + geom_bar(stat = 'identity', position = 'dodge')+
  scale_y_continuous(limits = c(-1,maximum+5),breaks =breaks)+
  labs(x = "",y = "Editing effeciency (%)")+
  annotate(geom="text", x=c(1:22), y=editrate+3, label=rates,
           color="black",size=3)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=10,colour='black'),axis.title=element_text(size=15,face="bold"),axis.text.x=element_text(angle=30,hjust=1,vjust=1))+
  guides(fill="none")
p
plotfile='./fig/on_target_editing_rate.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)
plotfile='./fig/on_target_editing_rate.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)



# venn plot ---------------------------------------------------------------
rm(list=ls())
library(ggplot2)
# install.packages('venn')
library(venn)
library(reshape2)
setwd('/data1/nyy/offtarget/v4')
mySNP<-read.table(file ='merge_26sample_filter_add_4method_overlap_venn.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mySNP)<-c('sample','SNP','method')
# unique(mySNP$sample)
# [1] "PC8"        "090202body" "cloningJR"  "PC7"        "PC5"        "PC6"        "cloningPF"  "PC4"       
# [9] "NC1"        "NC3"        "PC1"        "NC9"        "NC5"        "NC2"        "PC2"        "PC3"       
# [17] "NC6"        "cloningW"   "NC7"        "NC4"        "090202head" "NC8" 
samplelist<-c("090202body", "090202head",'NC1', 'NC2', 'NC3','NC4', 'NC5', 'NC6','NC7', 'NC8', 'NC9', 
              'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8',"cloningPF","cloningJR","cloningW")
# samplelist<-c('NC1')
for (i in samplelist){
  A1 = mySNP$SNP[mySNP$method == 'lofreq2' & mySNP$sample == i]
  B1 = mySNP$SNP[mySNP$method == 'platypus' & mySNP$sample == i]
  C1 = mySNP$SNP[mySNP$method == 'gatk' & mySNP$sample == i]
  D1 = mySNP$SNP[mySNP$method == 'strelka' & mySNP$sample == i]
  intersect1<-Reduce(intersect,list(A1,B1,C1,D1))
  print (intersect1)
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


# change count heatmap ----------------------------------------------------
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
options(scipen = 9)
mydata<-read.delim(file ='merge_26sample_filter_add_4method.txt',header = F, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('Individual','base_change','number')
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$Individual %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
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
plotfile='./fig/SNP.meancount.heatmap.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)
plotfile='./fig/SNP.meancount.heatmap.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)

# change  proportion heatmap ----------------------------------------------------
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
options(scipen = 9)

mydata<-read.delim(file ='merge_26sample_filter_add_4method_count_no_merge_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$Individual %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
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
plotfile='./fig/SNP.mean_proportion.heatmap.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)
plotfile='./fig/SNP.mean_proportion.heatmap.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)

############################################
rm(list=ls())
library(reshape2)
library(ggplot2)
library(dplyr)
options(scipen = 9)
mydata<-read.delim(file ='merge_26sample_filter_add_4method.txt',header = F, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('Individual','base_change','number')
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
mydata$Individual<-factor(mydata$Individual, levels=unique(mydata$Individual))

p<-ggplot(data = mydata, mapping = aes(x = factor(base_change), y = number,fill = Individual)) + geom_bar(stat = 'identity', position = 'dodge')+
  scale_y_continuous(limits = c(0,maximum),breaks =breaks)+
  labs(x = "Base change",y = "Number of SNVs")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=16),legend.title = element_blank())
p
plotfile='./fig/merge_26sample_filter_add_4method_count.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)
plotfile='./fig/merge_26sample_filter_add_4method_count.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 5)

#####
mydata<-read.delim(file ='merge_26sample_filter_add_4method_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
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
p
plotfile='./fig/merge_26sample_filter_add_4method_count_merge.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)
plotfile='./fig/merge_26sample_filter_add_4method_count_merge.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 5)


# #########boxplot --------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
options(scipen = 9)
mydata<-read.delim(file ='merge_26sample_filter_add_4method_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$Individual %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)

mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$base_change<-factor(mydata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
#mydata<-mydata[! mydata$Individual %in% c('cloningJR','cloningPF','cloningW','090202body','090202head','PC1','PC2','PC3','PC4'),]
#mydata<-mydata[! mydata$Individual %in% c('PC1','PC2','PC3','PC4'),]

breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]

p<-ggplot(mydata, aes(x=base_change, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  # stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum-1,hide.ns = TRUE,method.args = list(alternative = "less"))
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
# stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
#                    method = "wilcox.test",method.args =list(alternative = "greater"),size=5,label.y=maximum-1,hide.ns = TRUE)
p
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
  data1$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  data1$Individual %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  data1$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)

# shapiro.test(data1$number)
# 
# Shapiro-Wilk normality test
# 
# data:  data1$number
# W = 0.62671, p-value = 0.00000261


mydata<-rbind(data1,mydata)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$base_change<-factor(mydata$base_change, levels=c('ALL',"A>G","C>T","A>C","A>T","C>G","C>A"))
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]


p<-ggplot(mydata, aes(x=base_change, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('ALL','A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum-1,hide.ns = TRUE)
p
plotfile='./fig/base_change_barplot_PC_NC_v1_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/base_change_barplot_PC_NC_v1_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)


############only sum
sem <- function(x) sqrt(var(x)/length(x))
mydata<-mydata[mydata$base_change=='ALL',]
(ABE_mean<-mean(mydata$number[mydata$Status=='SCNT-ABE'])) ##124.8182
(ABE_sem<-sem(mydata$number[mydata$Status=='SCNT-ABE'])) ##84.45455
(SCNT_meam<-mean(mydata$number[mydata$Status=='SCNT']))##84.45455
(SCNT_sem<-sem(mydata$number[mydata$Status=='SCNT'])) ##6.932615

(ABE_sd<-sd(mydata$number[mydata$Status=='SCNT-ABE'])) ##107.4373
(SCNT_sd<-sd(mydata$number[mydata$Status=='SCNT']))   ##22.99288
(ABE_var<-var(mydata$number[mydata$Status=='SCNT-ABE'])) ##11542.76
(SCNT_var<-var(mydata$number[mydata$Status=='SCNT']))   ##528.6727


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


breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
mydata$Status
p<-ggplot(mydata, aes(x=Status, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  guides(color = "none", fill = "none")+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum-1,hide.ns = TRUE)

p
plotfile='./fig/4method_overlap_filter_sum_boxplot_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_filter_sum_boxplot_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

library("ggrepel")
breaks<-pretty(range(mydata$number), 8)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Status, y=number, fill=Status)) +
  geom_dotplot(binaxis='y', stackdir='center')+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="point")+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.title=element_text(size=20))+
  guides(color = "none", fill = "none")+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),
                     method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum+20,hide.ns = TRUE)

p1<-p+geom_text_repel(data=mydata, aes(label=Individual),col="black",alpha = 1)
p
p1

plotfile='./fig/4method_overlap_filter_sum_dotplot_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)
plotfile='./fig/4method_overlap_filter_sum_dotplot_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 5, height = 5)

plotfile='./fig/4method_overlap_filter_sum_dotplot_wilcoxlabel.pdf'
ggsave(plotfile, plot=p1, dpi = 600,width = 10, height = 6)
plotfile='./fig/4method_overlap_filter_sum_dotplot_wilcoxlabel.png'
ggsave(plotfile, plot=p1, dpi = 600,width = 10, height = 6)

# ##################################boxplot proportion --------------------
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_26sample_filter_add_4method_count_merge_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$Individual %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))

# mydata<-mydata[! mydata$Individual %in% c('PC1','PC2','PC3','PC4','cloningJR','cloningPF','cloningW','090202body','090202head'),]
# mydata<-mydata[! mydata$Individual %in% c('090202body'),]

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
plotfile='./fig/barplot_PC_NC_proportion_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/barplot_PC_NC_proportion_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)

######only SCNT-ABE
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_26sample_filter_add_4method_count_merge_proportion.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$Individual %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))

BEdata<-mydata[mydata$Status=='SCNT-ABE',]
BEdata$base_change<-factor(BEdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
# my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(BEdata, aes(x=base_change, y=proportion*100, color=base_change)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  # scale_y_continuous(limits = c(0,100),breaks = seq(0,100,10))+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/base_change_barplot_onlyPC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/base_change_barplot_onlyPC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

######only SCNT
CTdata<-mydata[mydata$Status=='SCNT',]
CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
# my_comparisons <- list(c('A>G', 'C>T'),c('C>T', "A>C"),c('C>T', "A>T"),c('C>T', "C>G"),c('C>T',"C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))

p<-ggplot(CTdata, aes(x=base_change, y=proportion*100, color=base_change)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  scale_y_continuous(limits = c(0,100),breaks = seq(0,100,20))+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/base_change_barplot_onlyNC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/base_change_barplot_onlyNC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

#############number
rm(list=ls())
library(ggplot2)
library(ggpubr)
options(scipen = 9)
mydata<-read.delim(file ='merge_26sample_filter_add_4method_count_merge.txt',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$Individual),]
mydata$Status<-case_when(
  mydata$Individual %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$Individual %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$Individual %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))

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
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  scale_x_discrete(labels=c('A>G/T>C','C>T/G>A','A>C/T>G','A>T/T>A','C>G/G>C','C>A/G>T'))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=16),legend.title = element_text(size=20))+
  guides(color=FALSE)+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/count_barplot_onlyPC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/count_barplot_onlyPC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

######NC
CTdata<-mydata[mydata$Status=='SCNT',]
CTdata$base_change<-factor(CTdata$base_change, levels=c("A>G","C>T","A>C","A>T","C>G","C>A"))
my_comparisons <- list(c('A>G', 'C>T'),c('A>G', "A>C"),c('A>G', "A>T"),c('A>G', "C>G"),c('A>G',"C>A"))
maximum<-max(CTdata$number)
maximum
breaks<-pretty(range(0,maximum*1.5), 8)
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
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test",label = "p.format",size=5,fontface= 'italic')
p
plotfile='./fig/count_barplot_onlyNC_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)
plotfile='./fig/count_barplot_onlyNC_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 8, height = 6)

#########gene_region -------------------------------------------------------------
####all_change type
rm(list=ls())
library(tidyverse)
mydata<-read.delim(file ='SNP_gene_region_count_all.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
mydata <- mydata %>%
  dplyr::group_by(sample) %>% 
  dplyr::mutate(proportion = number /sum(number))

mydata$Status<-case_when(
  mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$sample %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))

maximum<-max(mydata$number)
maximum
breaks<-pretty(range(0,maximum), 6)
maximum<- breaks[length(breaks)]
mydata$Region<-factor(mydata$Region, levels=c("UTR","Exonic","Intronic","Downstream","Upstream","Intergenic"))

p<-ggplot(mydata, aes(x=Region, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/SNP_gene_region_count_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/SNP_gene_region_count_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)

breaks<-pretty(range(0,1), 5)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Region, y=proportion, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=0.75,hide.ns = TRUE)
p
plotfile='./fig/SNP_gene_region_proportion_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/SNP_gene_region_proportion_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)


####################################only c('A>G','T>C'),]
rm(list=ls())
library(tidyverse)
mydata<-read.delim(file ='SNP_gene_region_count_AG.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
mydata <- mydata %>%
  dplyr::group_by(sample) %>% 
  dplyr::mutate(proportion = number /sum(number))

mydata$Status<-case_when(
  mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$sample %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$Region<-factor(mydata$Region, levels=c("UTR","Exonic","Intronic","Downstream","Upstream","Intergenic"))

maximum<-max(mydata$number)
maximum
breaks<-pretty(range(0,maximum), 6)
maximum<- breaks[length(breaks)]

p<-ggplot(mydata, aes(x=Region, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/SNP_gene_region_count_only_A2G_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/SNP_gene_region_count_only_A2G_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)

breaks<-pretty(range(0,1), 5)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Region, y=proportion, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/SNP_gene_region_proportion_only_A2G_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/SNP_gene_region_proportion_only_A2G_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)


####################################only c('C>T','G>A'),]
rm(list=ls())
library(tidyverse)
mydata<-read.delim(file ='SNP_gene_region_count_CT.txt',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
mydata<-mydata[order(mydata$sample),]
mydata <- mydata %>%
  dplyr::group_by(sample) %>% 
  dplyr::mutate(proportion = number /sum(number))

mydata$Status<-case_when(
  mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$sample %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$Status<-factor(mydata$Status, levels=c('SCNT-ABE', 'SCNT'))
mydata$Region<-factor(mydata$Region, levels=c("UTR","Exonic","Intronic","Downstream","Upstream","Intergenic"))

maximum<-max(mydata$number)
maximum
breaks<-pretty(range(0,maximum), 6)
maximum<- breaks[length(breaks)]

p<-ggplot(mydata, aes(x=Region, y=number, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Number of SNVs")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/SNP_gene_region_count_only_C2T_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/SNP_gene_region_count_only_C2T_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)

breaks<-pretty(range(0,1), 5)
maximum<- breaks[length(breaks)]
p<-ggplot(mydata, aes(x=Region, y=proportion, color=Status)) +
  geom_boxplot(position=position_dodge(0.8),outlier.colour = NA)+
  geom_jitter(position=position_dodge(0.8), size=0.8)+
  labs(x = "",y = "Percent of SNVs (%)")+
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour='black'),axis.ticks=element_line(size=1,colour='black'))+
  theme(axis.text=element_text(size=16,colour='black'),axis.text.x=element_text(angle=30,hjust=1,vjust=1),axis.title=element_text(size=20),legend.text=element_text(size=18),legend.title = element_blank())+
  stat_compare_means(aes(group = Status,label = paste0(sprintf("P=%.6s", ..p.format..),' ',..p.signif..)),method = "wilcox.test",size=5,fontface= 'italic',label.y=maximum,hide.ns = TRUE)
p
plotfile='./fig/SNP_gene_region_proportion_only_C2T_wilcox.pdf'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)
plotfile='./fig/SNP_gene_region_proportion_only_C2T_wilcox.png'
ggsave(plotfile, plot=p, dpi = 600,width = 10, height = 6)


############ggseqlogo AG
rm(list=ls())
# Load the required packages
setwd('/data1/nyy/offtarget/v4')
library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(stringi)
library(patchwork)
mydata<-read.delim(file ='all_edit_site_seqlog_AG.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq','changetype','sample')
mydata$Status<-case_when(
  mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$sample %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$sample1<-case_when(
  mydata$sample =="090202body"  ~  "SCNT10 (body)",
  mydata$sample =="090202head"  ~  "SCNT11 (head)",
  mydata$sample =="NC1"  ~  "SCNT1",
  mydata$sample =="NC2"  ~  "SCNT2",
  mydata$sample =="NC3"  ~  "SCNT3",
  mydata$sample =="NC4"  ~  "SCNT4",
  mydata$sample =="NC5"  ~  "SCNT5",
  mydata$sample =="NC6"  ~  "SCNT6",
  mydata$sample =="NC7"  ~  "SCNT7",
  mydata$sample =="NC8"  ~  "SCNT8",
  mydata$sample =="NC9"  ~  "SCNT9",
  mydata$sample =="cloningJR"  ~  "SCNT-ABE9 (Muscle)",
  mydata$sample =="cloningPF"  ~  "SCNT-ABE10 (Skin)",
  mydata$sample =="cloningW"  ~  "SCNT-ABE11 (Stomach)",
  mydata$sample =="PC1"  ~  "SCNT-ABE1",
  mydata$sample =="PC2"  ~  "SCNT-ABE2",
  mydata$sample =="PC3"  ~  "SCNT-ABE3",
  mydata$sample =="PC4"  ~  "SCNT-ABE4",
  mydata$sample =="PC5"  ~  "SCNT-ABE5",
  mydata$sample =="PC6"  ~  "SCNT-ABE6",
  mydata$sample =="PC7"  ~  "SCNT-ABE7",
  mydata$sample =="PC8"  ~  "SCNT-ABE8"
  )

samplelist<-c("SCNT1", "SCNT2", "SCNT3", "SCNT4", "SCNT5","SCNT6", "SCNT7", "SCNT8", "SCNT9",
              "SCNT10 (body)", "SCNT11 (head)", 
              "SCNT-ABE1", "SCNT-ABE2", "SCNT-ABE3", "SCNT-ABE4", "SCNT-ABE5", "SCNT-ABE6", "SCNT-ABE7", 
              "SCNT-ABE8", "SCNT-ABE9 (Muscle)", "SCNT-ABE10 (Skin)", "SCNT-ABE11 (Stomach)")

for (i in 1:22) {
  seq<-mydata$seq[mydata$sample1==samplelist[i]]
  p<-ggseqlogo(seq,method='p')+labs(title=samplelist[i])+theme(plot.title = element_text(hjust = 0.5,size=12))
  sname<-paste('s',i,sep = '')
  assign(sname, p)
  print(sname)
}

SCNTtext<-paste('s',1:11,sep = '')
SCNTABEtext<-paste('s',12:22,sep = '')

pSCNT<-eval(parse(text=paste(SCNTtext,collapse = '+')))
pSCNTABE<-eval(parse(text=paste(SCNTABEtext,collapse = '+')))

pall <- (pSCNT)/(pSCNTABE)   ####+plot_annotation(title = "", tag_levels = "A")
plotfile='./fig/all_edit_ggseqlogoAG.pdf'
ggsave(plotfile, plot=pall, dpi = 600,width = 8, height = 12)
plotfile='./fig/all_edit_ggseqlogoAG.png'
ggsave(plotfile, plot=pall, dpi = 600,width = 8, height = 12)


############ggseqlogo CT
rm(list=ls())
# Load the required packages
setwd('/data1/nyy/offtarget/v4')
library(ggplot2)
library(ggseqlogo)
library(ggpubr)
library(stringi)
library(patchwork)
mydata<-read.delim(file ='all_edit_site_seqlog_CT.txt',header = FALSE, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(mydata)<-c('pos','seq','changetype','sample')
mydata$Status<-case_when(
  mydata$sample %in% c('cloningJR','cloningPF','cloningW') ~ "SCNT-ABE",
  mydata$sample %in% c("090202body", "090202head",'NC1','NC2','NC3','NC4','NC5','NC6','NC7','NC8','NC9') ~ "SCNT",
  mydata$sample %in% c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8') ~ "SCNT-ABE"
)
mydata$sample1<-case_when(
  mydata$sample =="090202body"  ~  "SCNT10 (body)",
  mydata$sample =="090202head"  ~  "SCNT11 (head)",
  mydata$sample =="NC1"  ~  "SCNT1",
  mydata$sample =="NC2"  ~  "SCNT2",
  mydata$sample =="NC3"  ~  "SCNT3",
  mydata$sample =="NC4"  ~  "SCNT4",
  mydata$sample =="NC5"  ~  "SCNT5",
  mydata$sample =="NC6"  ~  "SCNT6",
  mydata$sample =="NC7"  ~  "SCNT7",
  mydata$sample =="NC8"  ~  "SCNT8",
  mydata$sample =="NC9"  ~  "SCNT9",
  mydata$sample =="cloningJR"  ~  "SCNT-ABE9 (Muscle)",
  mydata$sample =="cloningPF"  ~  "SCNT-ABE10 (Skin)",
  mydata$sample =="cloningW"  ~  "SCNT-ABE11 (Stomach)",
  mydata$sample =="PC1"  ~  "SCNT-ABE1",
  mydata$sample =="PC2"  ~  "SCNT-ABE2",
  mydata$sample =="PC3"  ~  "SCNT-ABE3",
  mydata$sample =="PC4"  ~  "SCNT-ABE4",
  mydata$sample =="PC5"  ~  "SCNT-ABE5",
  mydata$sample =="PC6"  ~  "SCNT-ABE6",
  mydata$sample =="PC7"  ~  "SCNT-ABE7",
  mydata$sample =="PC8"  ~  "SCNT-ABE8"
)

samplelist<-c("SCNT1", "SCNT2", "SCNT3", "SCNT4", "SCNT5","SCNT6", "SCNT7", "SCNT8", "SCNT9",
              "SCNT10 (body)", "SCNT11 (head)", 
              "SCNT-ABE1", "SCNT-ABE2", "SCNT-ABE3", "SCNT-ABE4", "SCNT-ABE5", "SCNT-ABE6", "SCNT-ABE7", 
              "SCNT-ABE8", "SCNT-ABE9 (Muscle)", "SCNT-ABE10 (Skin)", "SCNT-ABE11 (Stomach)")

for (i in 1:22) {
  seq<-mydata$seq[mydata$sample1==samplelist[i]]
  p<-ggseqlogo(seq,method='p')+labs(title=samplelist[i])+theme(plot.title = element_text(hjust = 0.5,size=12)) 
  sname<-paste('s',i,sep = '')
  assign(sname, p)
  print(sname)
}

SCNTtext<-paste('s',1:11,sep = '')
SCNTABEtext<-paste('s',12:22,sep = '')

pSCNT<-eval(parse(text=paste(SCNTtext,collapse = '+')))
pSCNTABE<-eval(parse(text=paste(SCNTABEtext,collapse = '+')))

pall <- (pSCNT)/(pSCNTABE)   ####+plot_annotation(title = "", tag_levels = "A")
plotfile='./fig/all_edit_ggseqlogoCT.pdf'
ggsave(plotfile, plot=pall, dpi = 600,width = 8, height = 12)
plotfile='./fig/all_edit_ggseqlogoCT.png'
ggsave(plotfile, plot=pall, dpi = 600,width = 8, height = 12)

