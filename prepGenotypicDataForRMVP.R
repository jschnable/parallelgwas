library(rMVP)
library(tidyverse)
args <- commandArgs(trailingOnly = FALSE)
infile <- str_remove(args[length(args)-1], fixed('-'))
outPrefix <- str_remove(args[length(args)], fixed('-'))
# load rmvp formatted data
MVP.Data(fileVCF=infile,
         fileKin=TRUE,
         filePC=TRUE,
         priority="memory",
         maxLine=10000,
         out=outPrefix
)
