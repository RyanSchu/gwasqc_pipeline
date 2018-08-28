library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--pihat", help="inbreeding coefficient to filter by")
parser$add_argument("--QCdir", help="directory where all the QC steps are written")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

##Following QCStep5b
ibd <- read.table(args$QCdir %&% "/QCStep5/QCStep5b/QCstep5b.genome", header = T)
pdf(args$QCdir %&% "/QCstats/IBD.pdf")
ggplot(data=ibd,aes(x=Z0,y=Z1))+geom_point(alpha=1/4)+theme_bw()
dev.off()
##pull duplicates
dups <- data.frame()
for(i in 1:dim(ibd)[1]){
  if(as.character(ibd$IID1[i]) == as.character(ibd$IID2[i])){
    dups <- rbind(dups,ibd[i,])
  }
}
dim(dups)
##Note and pull duplicates and missings
##In this example, there are neither duplicates nor missings
toExclude <- as.character(dups$IID1)
a <- as.character(ibd$IID1) %in% toExclude
others <- ibd[a==FALSE,]
dim(others)

pdf(args$QCdir %&% "/QCstats/ibd_pi_hat.pdf")
hist(others$PI_HAT)
dev.off()

#sortOthers <- others[order(others$PI_HAT, decreasing = TRUE),]
##Unexpected duplicates:
#filter(others,PI_HAT>=0.2)
write.table(filter(others,PI_HAT>=args$pihat), args$QCdir %&% "/QCstats/related.to.remove.txt", quote = FALSE, row.names = FALSE)
