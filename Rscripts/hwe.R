library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--hwe", help="full path to the plink .hwe file to be analyzed")
parser$add_argument("-t", "--tag", help="name tag for identification")
parser$add_argument("-o", "--outputdir", help="directory where you would like to output your plots")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

hwe <- read.table(args$hwe, header = T)
pdfname<-paste(args$outputdir,"/hwestats", args$tag, ".pdf", sep="")
pdf(pdfname)
hist(hwe$P)
dev.off()

hwestats<-paste(args$outputdir,"/hwestats", args$tag, ".txt", sep="")
write("Min. 1st Qu. Median Mean 3rd Qu. Max.", hwestats)
write(summary(hwe$P),hwestats, append=T)
write("P<1e-06   P>1e-06",hwestats,append=T)
write(table(hwe$P<1e-06), hwestats, append=T)
write(table(hwe$P<1e-06)/sum(table(hwe$P<1e-06)), hwestats, append=T)

