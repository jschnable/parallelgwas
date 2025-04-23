library(readr)
library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
infile <- str_remove(args[length(args)], fixed('-'))

processLD <- function(ld_file)
{
  ld <- read.table(ld_file, header = TRUE)

  ld <- ld %>%
    rowwise() %>%
    mutate(snp1_bp = min(BP_A, BP_B),
           snp2_bp = max(BP_A, BP_B),
           snp1 = case_when(BP_A < BP_B ~ SNP_A, .default = SNP_B),
           snp2 = case_when(BP_A < BP_B ~ SNP_B, .default = SNP_A)) %>%
    select(snp1, snp1_bp, snp2, snp2_bp)
  return(ld)
}

write_csv(processLD(infile), str_replace(infile, '.ld', '.csv'), quote = 'needed')
