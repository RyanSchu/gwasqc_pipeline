#read in diff chr results
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--column","-c", help="Specifies which column contains pvalues for multi-column file", type="integer")
parser$add_argument("--input","-i",help="file containing pvalues. Can be fed one column file as pvalues, or can specify column number for multi-column file with '--column'. If neither of these is the case then it assumes it is plink format by default.")
parser$add_argument("--limit", help="Limits output to the top proportion of data you wish to graph, default is no limit. Accepts decimal input 0 to 1.", type="double")
parser$add_argument("--noheader", help="Specifies that input file does not contain a header. Default assumes file has header", action='store_true')
parser$add_argument("--out", "-o", help="Use as in plink - output prefix optionally including file path, but not including file type.")
parser$add_argument("--range", help="Specifies the total number of pvalues within the parent set of input. For use when input pvalues is a subset of the total.", type="integer")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

##three situations - standard plink output, only one column, or specify column number
##check if theres a header or not

read_input<- function(file=args$input, noheader=args$noheader, limit=Inf) {
  if (noheader == T && grepl(file,".gz$") == T){
    print("reading input")
    input<-as.data.frame(fread('zcat ' %&% file, header = F, showProgress = T))
    return(input)
  } else if (noheader == F && grepl(file,".gz$") == T){
    print("reading input")
    input<-as.data.frame(fread('zcat ' %&% file, header = T, showProgress = T))
    return(input)
  } else if (noheader == F && grepl(file,".gz$") == F){
    print("reading input")
    input<-as.data.frame(fread(file, header = T, showProgress = T))
    return(input)
  } else if (noheader == T && grepl(file,".gz$") == F){
    print("reading input")
    input<-as.data.frame(fread(file, header = F, showProgress = T))
    return(input)
  }
}

clean_pvals<-function(df){ #grabs the pvalues, -log transforms, sorts them, drops NAs
  if (ncol(df) == 1){ #if only one col, that is just assumed as a list of pvals
    print("Input data has only one column, treating it as pvalues")
    logp<- -log10(df)
    pval_col<- sort(logp, decreasing = T, na.last = NA) %>% as.data.frame()
    return(pval_col)
  }
  else if(!is.null(args$column)){ #have you specified which col number is pvalues
    print("Treating i-th column as pvalues where i = " %&% args$column)
    logp<- -log10(df[,args$column])
    pval_col<-sort(logp, decreasing = T, na.last = NA) %>% as.data.frame()
    return(pval_col)
  }
  else { #default assumes it is plink format
    colnames(df)<-c("CHR","SNP","BP","A1","TEST","NMISS","OR","STAT","P")
    pval_col<-filter(df, Test == "ADD") %>% select(P) %>% -log10() %>% sort(decreasing = T, na.last = NA) %>% as.data.frame()
    return(pval_col)
  }
}

create_qq_input<-function(df, limit=nrow(df),range=nrow(df)){ #takes in one column df of pvalues - assumes no NAs, can also limit to the top n hits. Default of this will be just all pvalues
  print("Size of parental set is presumed to be " %&% range)
  print("Provided number of pvalues is " %&% nrow(df))
  print("limiting output to top " %&% limit %&% " rows")
  count <- nrow(df) #count number of pvals that are not na
  ExpP <- -log10((1:count)/(range+1)) %>% as.data.frame() #generate a list of expected pvalues equal in length to the obs pvalues that are not na
  qqvals<-cbind.data.frame(df[1:limit,],ExpP[1:limit,])
  colnames(qqvals)<-c("Observed","Expected")
  return(qqvals)
}

pvalue_data<-read_input(args$input, args$noheader)
cleanded_pvalues<-clean_pvals(pvalue_data)
#set limit
if (!is.null(args$limit)){
  limit = floor(args$limit * nrow(cleanded_pvalues))
} else {
  limit = nrow(cleanded_pvalues)
}
if (!is.null(args$range)){
  range = args$range
} else {
  range = nrow(cleanded_pvalues)
  
}
qq_input<-create_qq_input(cleanded_pvalues, limit, range)

qq1 <- ggplot(qq_input,aes(x=Expected,y=Observed)) + 
  geom_point(shape=1) + 
  #coord_cartesian(xlim=c(-0.05,8.05),ylim=c(-0.05,30.05)) + 
  theme_bw(12) + 
  scale_colour_manual(values=c(rgb(163,0,66,maxColorValue = 255),"dark gray")) + 
  theme(legend.position = c(0.01,0.99),legend.justification = c(0,1)) + 
  geom_abline(slope=1,intercept = 0) + 
  xlab(expression(Expected~~-log[10](italic(p)))) + 
  ylab(expression(Observed~~-log[10](italic(p)))) +
  expand_limits(x = 0, y = 0)


png(filename = args$out %&% '.png', width=480, height=480)
qq1
dev.off()

pdf(file = args$out %&% '.pdf',width = 3.75, height = 3.75)
qq1
dev.off()

