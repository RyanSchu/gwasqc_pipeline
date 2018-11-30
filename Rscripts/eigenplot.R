library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--val",help="full path to eigenvalue file")
parser$add_argument("--vec",help="full path to eigenvector file")
parser$add_argument("-o", "--outputdir", help="directory where you would like to output your plots")
args <- parser$parse_args()
"%&%" = function(a,b) paste (a,b,sep="")

pcaplots <- args$outputdir %&% "/unmerged_pca_plots.pdf"
pcs <- read.table(args$vec,header=T)
eval <- scan(args$val)[1:10]
skree<-round(eval/sum(eval),3)
skree<-cbind.data.frame(skree,c(1,2,3,4,5,6,7,8,9,10))
colnames(skree)<-c("percent_var", "PC")

pdf(pcaplots)

ggplot(data=skree, aes(x=PC, y=percent_var)) + geom_point() + geom_line() + scale_x_continuous(breaks = 1:10) + ggtitle("Proportion of variance explained")

#PCA Plot 1 (PC1 vs PC2)
ggplot() + geom_point(data=pcs,aes(x=PC1,y=PC2)) + theme_bw() + scale_colour_brewer(palette="Set1") + ggtitle("PC1 vs PC2")

#PCA Plot 2 (PC1 vs PC3)
ggplot() + geom_point(data=pcs,aes(x=PC1,y=PC3)) + theme_bw() + scale_colour_brewer(palette="Set1") + ggtitle("PC1 vs PC3")

#PCA Plot 1 (PC2 vs PC3)
ggplot() + geom_point(data=pcs,aes(x=PC2,y=PC3)) + theme_bw() + scale_colour_brewer(palette="Set1") + ggtitle("PC2 vs PC3")

dev.off() 
