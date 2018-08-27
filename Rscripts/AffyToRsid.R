library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--bim", help="full path to the bim file you would like to translate")
parser$add_argument("--csv",help="full path to file containing the Affy SNPs with corresponding rsids")
parser$add_argument("--QCdir", help="directory where all the QC steps are written")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

csv<-read.table(args$csv, header = T, skip = 18, sep = ',')
bim<-read.table(args$bim, sep = '\t')
bim<- left_join(bim, csv, by = c('V2' = 'Probe.Set.ID'))
newbim<- select(bim, c('V1', 'dbSNP.RS.ID', 'V3', 'V4', 'V5', 'V6' ))
newbim<- newbim[!grepl("---", newbim$dbSNP.RS.ID),]
newbim<- na.omit(newbim)

write.table(x= newbim, file = args$bim, quote = F, row.names = F,col.names = F)
