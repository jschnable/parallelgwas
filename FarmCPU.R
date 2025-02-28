local({
    r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos = r)
})
if (!require(rMVP)){ install.packages('rMVP') }

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

for(i in 2:ncol(pheno))
{
  imMVP<-MVP(phe = pheno[,c(1,i)], geno = genotype, map = map, K=Kinship, CV.FarmCPU=Covariates_PC,
               nPC.FarmCPU = 3, maxLoop = 10, method = "FarmCPU", p.threshold = (0.05/EFFECTIVE_MARKERS), 
               threshold = (0.05/effective_ratio),
               file.output = 'pmap.signal')
}

