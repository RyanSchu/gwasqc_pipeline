library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-g", "--genome", help="Full path to the .genome file you wish to use")
parser$add_argument("-o", "--outputdir", help="directory where you would like to output your plots")
#parser$add_argument("-p", "--pihat", help="inbreeding coefficient to filter by")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

#Following QCStep5b
ibd <- read.table(args$genome, header = T)
png(args$outputdir %&% "/IBD.png")
ggplot(data=ibd,aes(x=Z0,y=Z1))+geom_point(alpha=1/4)+theme_bw() + xlim(0,1) + ylim(0,1)
dev.off()
#pull duplicates
dups <- data.frame()
for(i in 1:dim(ibd)[1]){
  if(as.character(ibd$IID1[i]) == as.character(ibd$IID2[i])){
    dups <- rbind(dups,ibd[i,])
  }
}
dim(dups)
#Note and pull duplicates and missings
#In this example, there are neither duplicates nor missings
toExclude <- as.character(dups$IID1)
a <- as.character(ibd$IID1) %in% toExclude
others <- ibd[a==FALSE,]
dim(others)

pdf(args$outputdir %&% "/ibd_pi_hat.pdf")
hist(others$PI_HAT)
dev.off()

#sortOthers <- arrange(others, PI_HAT, )
#Unexpected duplicates:
#filter(others,PI_HAT>=0.2)
write.table(filter(others,PI_HAT>=0.25), args$outputdir %&% "/related.and.duplicates.txt", quote = FALSE, row.names = FALSE)
