#devtools::install_github("argparse", "trevorld")
library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-t", "--threshold", help="whatever threshold you filtered by", type="double", default=0.01 )
parser$add_argument("--QCdir", help="directory where all the QC steps are written")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")
percentage<-100*(1-args$threshold)

#read in files
lmiss <- read.table(args$QCdir %&% "/missingness_hwe_steps/01initial_missingness.lmiss", header=T)
imiss <- read.table(args$QCdir %&% "/missingness_hwe_steps/01initial_missingness.imiss", header=T)
newimiss <- read.table(args$QCdir %&% "/missingness_hwe_steps/03missingness_validation.imiss", header = T)
newlmiss <- read.table(args$QCdir %&% "/missingness_hwe_steps/03missingness_validation.lmiss", header = T)

#create histogram's of missingness before and after filtering
pdf(args$QCdir %&% "/plots_stats/callRateDistributions.pdf")
hist(lmiss$F_MISS, main="Call rate distribution for SNPs before filtering")
hist(newlmiss$F_MISS, main="distribution of SNPs after removing SNPs with call rate < " %&% percentage)
hist(imiss$F_MISS, main="distribution of individuals before filtering")
hist(newimiss$F_MISS, main="distribution of individuals after SNPs filtering")
dev.off()

#stats
statsfile<- args$QCdir %&% "/plots_stats/missingness.txt"
#Initial SNP count, SNP count after filtering, individual count after filtering
init<- "Initial number of SNPs is " %&% dim(lmiss)[1]
final<- paste("Number of SNPs after filtering with out call rates < ", percentage, " is ", dim(newlmiss)[1], sep="")
indiv<- "Number of individuals after filtering is " %&% dim(imiss)[1]
write(init, statsfile, append=T)
write(final, statsfile, append=T)
write(indiv, statsfile, append=T)




