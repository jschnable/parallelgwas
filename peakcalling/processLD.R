library(readr)
library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
infile <- str_remove(args[length(args)], fixed('-'))

processLD <- function(ld_file)
{
  ld <- read.table(ld_file, header = TRUE)
  
  ld <- ld %>% 
    mutate(ordered = BP_A < BP_B)
  
  ld_ordered <- ld %>%
    filter(ordered) %>% 
    rename(snp1_bp = BP_A, 
           snp1 = SNP_A, 
           snp2_bp = BP_B, 
           snp2 = SNP_B)
  
  ld_reversed <- ld %>% 
    filter(!ordered) %>%
    rename(snp1_bp = BP_B, 
           snp1 = SNP_B,
           snp2_bp = BP_A, 
           snp2 = SNP_A)

  ld <- bind_rows(ld_ordered, ld_reversed)
  return(ld)
}

write_csv(processLD(infile), str_replace(infile, '.ld', '.csv'), quote = 'needed')
