library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--QCdir", help="directory where all the QC steps are written")
parser$add_argument("-t", "--threshold", help="call rate filtering threshold. Keeps SNPs with call rates > 1-threshold.", type="double", default=0.01 )
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

pdf(args$QCdir %&% "/QCplots/callRateDistributions.pdf")
percentage<-100*(1-args$threshold)

##Following QCStep1
##displays the distribution of proportion of SNPs (F_MISS) missing in the sample
lmiss <- read.table(args$QCdir %&% "/QCStep1/QCStep1.lmiss", header=T)
hist(lmiss$F_MISS, main="Call rate distribution for SNPs before filtering")
##Beginning SNP count
dim(lmiss)[1]
##SNPs with call rates > 99%
table(lmiss$F_MISS<0.01)
##percent SNPs with call rates > 99%
table(lmiss$F_MISS<0.01)/sum(table(lmiss$F_MISS<0.01))

##Following QCStep3
imiss <- read.table(args$QCdir %&% "/QCStep3/QCStep3.imiss", header = T)
hist(imiss$F_MISS, main="distributions for individuals after removing SNPs call rate < " %&% percentage)
##Now plot SNPs with > 99% call rates
newlmiss <- read.table(args$QCdir %&% "/QCStep3/QCStep3.lmiss", header = T)
hist(newlmiss$F_MISS, main="distributions for SNPs after removing SNPs call rate < " %&% percentage)

##SNP and individual count after rm low-call SNPs
dim(newlmiss)[1]
dim(imiss)[1]

dev.off()
