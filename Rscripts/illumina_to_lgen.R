library(dplyr)
library(tidyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--illumina", help="file path of the illumina FinalReport.txt")
parser$add_argument("-o", "--out", help="Use as in plink. Out is output prefix including full path and no file type extensions.")
parser$add_argument("--skip", help="number of lines that preceed the header in input file. Default is 9", type="integer", default=9)
args <- parser$parse_args()

#Process the final report.txt. Assumes that there are 9 header lines
Finalreport<-as.data.frame(read.table(file=args$illumina, sep='\t', skip = args$skip, header = T))
Finalreport<-filter(Finalreport, Allele1...Forward != 'I')
Finalreport["empty"]="0"
Finalreport['fid']<-Finalreport$Sample.ID

#create map df
map<-select(Finalreport, Chr, SNP.Name, empty, Position)
map<-map[!duplicated(map),]
map<-map[complete.cases(map),]

#create lgen df
lgen<-select(Finalreport, fid, Sample.ID, SNP.Name, Allele1...Forward, Allele2...Forward)
lgen<-lgen[!duplicated(lgen),]
lgen<-lgen[complete.cases(lgen),]
lgen<-filter(lgen, Allele1...Forward != "-" & Allele2...Forward != "-")

#write to file
write.table(map, file = paste(args$out, ".map",sep=""), sep = "\t", col.names = F, row.names = F, quote = F)
write.table(lgen, file = paste(args$out,".lgen",sep=""), sep = "\t", col.names = F, row.names = F, quote = F)
