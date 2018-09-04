library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--het", help="full path to the het file")
parser$add_argument("--tag", "-t", help="descriptive tag for file naming")
parser$add_argument("-o", "--outputdir", help="directory where you would like to output your plots")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

hetoutfile<-paste(args$outputdir,"/Hetoutliers", args$tag, ".txt", sep="")
hetpdf<-paste(args$outputdir,"/Heterozygosity", args$tag, ".pdf", sep="")
hetstats<-paste(args$outputdir,"/Hetstats", args$tag, ".txt", sep="")

HET <- read.table(args$het, header = T, as.is = T)
H = (HET$N.NM.-HET$O.HOM.)/HET$N.NM.
oldpar = par(mfrow=c(1,2))

pdf(hetpdf)
hist(H,50, main="H estimate")
hist(HET$F,50, main="Heterozygosity estimates")
abline(v=mean(HET$F)+6*sd(HET$F),col="red")
abline(v=mean(HET$F)-6*sd(HET$F),col="red")
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
write.table(allexclude2, file = hetoutfile, quote = F, col.names = F, row.names = F)

