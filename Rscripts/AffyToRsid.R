library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--bim", help="full path to the bim file you would like to translate")
parser$add_argument("-o", "--outputdir", help="directory where you would like to output your plots")
parser$add_argument("--stats",help="full path to the stats file containing the rsIDs")
args <- parser$parse_args()

"%&%" = function(a,b) paste (a,b,sep="")

#csv<-read.table(args$csv, header = T, skip = 18, sep = ',')
#bim<-read.table(args$bim, sep = '\t')
#bim<- left_join(bim, csv, by = c('V2' = 'Probe.Set.ID'))
#newbim<- select(bim, c('V1', 'dbSNP.RS.ID', 'V3', 'V4', 'V5', 'V6' ))
#newbim<- newbim[!grepl("---", newbim$dbSNP.RS.ID),]
#newbim<- na.omit(newbim)

#write.table(x= newbim, file = args$bim, quote = F, row.names = F,col.names = F)


#Dependent on the genotyping platoform used (i.e. Affymatrix, Illumina, etc.),
#the SNP identifiers are recorded differently (rsID, SNP_A-#, AFFY-SNP-#).
#The hapmap individuals are recorded by rsID, and they share different positions from SNP_A-#.
#This means that we have to change the MGS data identifiers to rsID.

TotalSNPs<-read.table(args$stats, sep="\t", header=T)
#The summary statistics from the dbGaP data containing information that contains an rsID and position for each SNP_A-#
SelectSNPs<-dplyr::select(TotalSNPs,"MarkerAccession","ChrID", "ChrPosition","SubmittedSNPID")
#Isolates the positions we need to merge with a .bim file
bim<-read.table(args$bim, header = F)
#Reading in the last .bim from QC

mergedbim <- left_join(bim, SelectSNPs, by = c("V1" = "ChrID", "V4" = "ChrPosition"))
mergedbim2 <- mutate(mergedbim,snp=ifelse(is.na(MarkerAccession), as.character(V2), as.character(MarkerAccession)))
newbim <- dplyr::select(mergedbim2,V1,snp,V3,V4,V5,V6)
filetest2<-newbim[!duplicated.data.frame(newbim),]

write.table(filetest2,file = args$outputdir %&% "/00rsID_format.bim",quote=F, sep="\t",row.names=F,col.names=F)
