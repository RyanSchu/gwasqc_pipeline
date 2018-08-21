library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-r","--relatedness",help="Threshold of relatedness used for filtering")
parser$add_argument("--QCdir", help="directory where all the QC steps are written")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

hetstats<-args$QCdir %&% "/QCstats/Hetstats.txt"
HET <- read.table(args$QCdir %&% "/QCStep5/QCStep5c/QCStep5c.het", header = T, as.is = T)
newHET <- read.table(args$QCdir %&% "/QCStep5/QCStep5e/QCStep5e.het", header = T, as.is = T)
H = (HET$N.NM.-HET$O.HOM.)/HET$N.NM.
oldpar = par(mfrow=c(1,2))

pdf(args$QCdir %&% "/QCstats/Heterozygosity.pdf")
hist(H,50, main="H estimate before filtering")
hist(HET$F,50, main="Heterozygosity estimates prior to filtering")
abline(v=mean(HET$F)+6*sd(HET$F),col="red")
abline(v=mean(HET$F)-6*sd(HET$F),col="red")
hist(newHET$F,50, main="Heterozygosity estimates after filtering by relatedness threshold of " %&% args$relatedness)
abline(v=mean(newHET$F)+6*sd(newHET$F),col="red")
abline(v=mean(newHET$F)-6*sd(newHET$F),col="red")
dev.off()

write("Min. 1st Qu. Median Mean 3rd Qu. Max. of HET$F", hetstats)
write(summary(HET$F),hetstats, append=T)

par(oldpar)

sortHET <- HET[order(HET$F),]
outliers <- data.frame()

for(i in 1:length(sortHET$F)){
  if(sortHET[i,6] > (mean(sortHET$F)+3*sd(sortHET$F))){
    outliers <- rbind(outliers, sortHET[i,])
  }
  if(sortHET[i,6] < (mean(sortHET$F)-3*sd(sortHET$F))){
    outliers <- rbind(outliers, sortHET[i,])
  }
}
hetoutliers <- select(outliers, FID, IID)
dim(hetoutliers)
allexclude2 <- hetoutliers
write.table(allexclude2, file = args$QCdir %&% "/QCStep5/QCStep5c/QCStep5c.txt", quote = F, col.names = F, row.names = F)


newsortHET <- newHET[order(newHET$F),]
newoutliers <- data.frame()

for(i in 1:length(newsortHET$F)){
  if(newsortHET[i,6] > (mean(newsortHET$F)+3*sd(newsortHET$F))){
    newoutliers <- rbind(newoutliers, newsortHET[i,])
  }
  if(newsortHET[i,6] < (mean(newsortHET$F)-3*sd(newsortHET$F))){
   newoutliers <- rbind(newoutliers, newsortHET[i,])
  }
}
newhetoutliers <- select(outliers, FID, IID)
dim(newhetoutliers)
newallexclude2 <- newhetoutliers
write.table(newallexclude2, file = args$QCdir %&% "/QCStep5/QCStep5e/QCStep5e.txt", quote = F, col.names = F, row.names = F)
count <- paste(dim(newhetoutliers)[1], "to remove after filtering")
write(count, hetstats, append=T)
