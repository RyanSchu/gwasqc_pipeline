library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--QCdir", help="directory where all the QC steps are written")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")


hwe <- read.table(args$QCdir %&% "/QCStep4/QCStep4.hwe", header = T)
summary(hwe$P)
pdf(args$QCdir %&% "/QCstats/hwestats.pdf")
hist(hwe$P)
dev.off()

table(hwe$P<1e-06)
table(hwe$P<1e-06)/sum(table(hwe$P<1e-06))
``
