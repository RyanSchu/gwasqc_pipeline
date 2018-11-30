library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--hapmapdir", help="directory where all the hapmap files are written")
parser$add_argument("--val",help="full path to eigenvalue file")
parser$add_argument("--vec",help="full path to eigenvector file")
parser$add_argument("--fam", help="full path to the fam file youd like to use")
parser$add_argument("-o", "--outputdir", help="directory where you would like to output your plots")
parser$add_argument("--pop", help="full path to the ")
args <- parser$parse_args()
"%&%" = function(a,b) paste (a,b,sep="")

pcaplots <- args$outputdir %&% "/merged_pca_plots.pdf"

if (!is.null(args$pop)){
  hapmappopinfo <- read.table(args$pop) %>% select (V1,V3)
} else if (grepl("19",args$hapmapdir, fixed=TRUE)) {
  hapmappopinfo <- read.table(args$hapmapdir %&% "/pop_HM3_hg19_forPCA.txt") %>% select (V1,V3)
} else if (grepl( "18",args$hapmapdir, fixed=TRUE)) {
  hapmappopinfo <- read.table(args$hapmapdir %&% "/pop_HM3_hg18_forPCA.txt") %>% select (V1,V3)
}
colnames(hapmappopinfo) <- c("pop","IID")
fam <- read.table(args$fam) %>% select (V1,V2)
colnames(fam) <- c("FID","IID")
popinfo <- left_join(fam,hapmappopinfo,by="IID")
popinfo <- mutate(popinfo, pop=ifelse(is.na(pop),'GWAS', as.character(pop)))
table(popinfo$pop)
pcs <- read.table(args$vec,header=T)
pcdf <- data.frame(popinfo, pcs[,3:ncol(pcs)])
gwas <- filter(pcdf,pop=='GWAS')
hm3 <- filter(pcdf, grepl('NA',IID))
eval <- scan(args$val)[1:10]
skree<-round(eval/sum(eval),3)
skree<-cbind.data.frame(skree,c(1,2,3,4,5,6,7,8,9,10))
colnames(skree)<-c("percent_var", "PC")

pdf(pcaplots)

ggplot(data=skree, aes(x=PC, y=percent_var)) + geom_point() + geom_line() + scale_x_continuous(breaks = 1:10) + ggtitle("Proportion of variance explained")

#PCA Plot 1 (PC1 vs PC2)
ggplot() + geom_point(data=pcdf,aes(x=PC1,y=PC2,col=pop,shape=pop)) + theme_bw() + scale_colour_brewer(palette="Set1") + ggtitle("PC1 vs PC2")

#PCA Plot 2 (PC1 vs PC3)
ggplot() + geom_point(data=pcdf,aes(x=PC1,y=PC3,col=pop,shape=pop)) + theme_bw() + scale_colour_brewer(palette="Set1") + ggtitle("PC1 vs PC3")

#PCA Plot 1 (PC2 vs PC3)
ggplot() + geom_point(data=pcdf,aes(x=PC2,y=PC3,col=pop,shape=pop)) + theme_bw() + scale_colour_brewer(palette="Set1") + ggtitle("PC2 vs PC3")

#PCA with HAPMAP populations
yri <- filter(pcdf,pop=='YRI')
uPC1 <- mean(yri$PC1) + 5*sd(yri$PC1)
lPC1 <- mean(yri$PC1) - 5*sd(yri$PC1)
uPC2 <- mean(yri$PC2) + 5*sd(yri$PC2)
lPC2 <- mean(yri$PC2) - 5*sd(yri$PC2)
ggplot() + geom_point(data=gwas,aes(x=PC1,y=PC2,col=pop,shape=pop))+geom_point(data=hm3,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw() +geom_vline(xintercept=c(uPC1,lPC1)) +geom_hline(yintercept=c(uPC2,lPC2)) + ggtitle("Assuming homogeneous, non-admixed")



inclusion <- gwas[gwas$PC1 >= lPC1,]
inclusion <- inclusion[inclusion$PC1 <= uPC1,]
inclusion <- inclusion[inclusion$PC2 >= lPC2,]
inclusion <- inclusion[inclusion$PC2 <= uPC2,]
samples <- inclusion[,1:2]
table(inclusion$pop)

dim(samples)[1]
dim(gwas)[1]-dim(samples)[1]

ggplot() + geom_point(data=gwas,aes(x=PC1,y=PC2,col=gwas$IID %in% samples$IID,shape=gwas$IID %in% samples$IID))+geom_point(data=hm3,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw() + ggtitle("Assuming homogeneous, non-admixed")

dev.off() 
#write.table(samples, args$QCdir %&% "/PCA/GWAS_PCA.txt",quote=F,row.names=F,col.names=F)

#afrpcs <- read.table("/home/angela/px_yri_chol/QC/QCStep6/QCStep6e/QCStep6e.evec",skip=1)
#afrcdf <- afrpcs %>% rename(PC1=V2,PC2=V3,PC3=V4,PC4=V5,PC5=V6,PC6=V7,PC7=V8,PC8=V9,PC9=V10,PC10=V11) %>% mutate(pop=ifelse(grepl("TC",V1),"GWAS","GWAS"))
#eval <- scan("/home/angela/px_yri_chol/QC/QCStep6/QCStep6e/QCStep6e.eval")[1:10]
#round(eval/sum(eval),3)
