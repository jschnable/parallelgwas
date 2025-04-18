library(rMVP)
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
pheno_file <- str_remove(args[length(args)-3], fixed('-'))
geno_stem <- str_remove(args[length(args)-2], fixed('-'))
TOTAL_MARKERS <- as.numeric(str_remove(args[length(args)-1], fixed('-')))
EFFECTIVE_MARKERS <- as.numeric(str_remove(args[length(args)], fixed('-'))) 

effective_ratio <- EFFECTIVE_MARKERS/TOTAL_MARKERS
pheno <- read.csv(pheno_file)
genotype <- attach.big.matrix(str_c(geno_stem, ".geno.desc"))
map <- read.table(str_c(geno_stem, ".geno.map"), header = TRUE)
Kinship <- attach.big.matrix(str_c(geno_stem, ".kin.desc"))
Covariates_PC <- bigmemory::as.matrix(attach.big.matrix(str_c(geno_stem, ".pc.desc")))

Sys.setenv(OMP_NUM_THREADS=1)

# Load individuals file
ind <- read_delim(str_c(geno_stem, ".geno.ind"), delim = '\t', col_names = FALSE)

pheno <- left_join(ind, pheno, join_by(X1==genotype)) %>% rename(genotype = X1)

for(i in 2:ncol(pheno))
{
  imMVP<-MVP(phe = pheno[,c(1,i)], geno = genotype, map = map, K=Kinship, CV.MLM=Covariates_PC,
               nPC.MLM = 3, maxLine=10000, method = "MLM", method.bin="static", p.threshold = (0.05/EFFECTIVE_MARKERS), 
               threshold = (0.05/effective_ratio),
               file.output = 'pmap.signal', vc.method="BRENT", ncpus=15)
}

