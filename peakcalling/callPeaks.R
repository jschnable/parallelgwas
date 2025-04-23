library(readr)
library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
ldfile <- str_remove(args[length(args)-4], fixed('-'))
signalfile <- str_remove(args[length(args)-3], fixed('-'))
chrom <- as.numeric(str_remove(args[length(args)-2], fixed('-')))
traitlist <- str_remove(args[length(args)-1], fixed('-'))
tag <- str_remove(args[length(args)], fixed('-'))

callLocalPeaks <- function(signals, ld, chr, trait, min_snps_per_peak=3)
{
  n_snps <- length(signals$SNP)
  if(n_snps < min_snps_per_peak)
  {
    return(NULL)
  }
  peak_df <- tibble()
  signals <- signals %>% 
    arrange(POS)
  ld <- ld %>% 
    filter((snp1 %in% unique(signals$SNP) & (snp2 %in% signals$SNP))) %>% 
             arrange(snp1_bp, snp2_bp)
  # set initial values
  cur_peak_start <- signals$POS[1]
  cur_peak_end <- cur_peak_start
  cur_peak_n_snps <- 1
  cur_peak_min_pval <- signals$pval[1]
  cur_peak_top_snp <- signals$SNP[1]
  i <- 1
  
  while(i <= n_snps)
  {
    if(n_snps - i + cur_peak_n_snps < min_snps_per_peak)
    {
      if(nrow(peak_df > 0))
      {
        return(peak_df)
      }
      else
      {
        print(str_c('No peaks found on chromosome ', chr, ' for ', trait))
        return()
      }
    }
    
    cur_ld <- filter(ld, snp1==signals$SNP[i])
    n_snps_in_ld <- nrow(cur_ld)
    # nothing in LD, 
    if(n_snps_in_ld==0)
    {
      #and we don't have enough snps before this to call it a peak: update current peak values to next snp
      if(cur_peak_n_snps < min_snps_per_peak)
      {
        cur_peak_start <- signals$POS[i+1]
        cur_peak_end <- cur_peak_start
        cur_peak_min_pval <- signals$pval[i+1]
        cur_peak_top_snp <- signals$SNP[i+1]
        # and reset number of snps in peak
        cur_peak_n_snps <- 1
        i <- i + 1
      }
      # end  peak
      else
      {
        if(signals$pval[i] < cur_peak_min_pval)
        {
          cur_peak <- tibble_row(start = cur_peak_start, 
                                 end = signals$POS[i], 
                                 min_pval = signals$pval[i], 
                                 top_snp = signals$SNP[i],
                                 n_snps = cur_peak_n_snps,
                                 chrom = chr,
                                 phenotype = trait)
          peak_df <- bind_rows(peak_df, cur_peak)
        }
        else
        {
          cur_peak <- tibble_row(start = cur_peak_start, 
                                 end = signals$POS[i], 
                                 min_pval = cur_peak_min_pval, 
                                 top_snp = cur_peak_top_snp,
                                 n_snps = cur_peak_n_snps,
				 chrom = chr,
           			 phenotype = trait)
          peak_df <- bind_rows(peak_df, cur_peak)
        }
        # update current values to next snp
        cur_peak_start <- signals$POS[i+1]
        cur_peak_end <- cur_peak_start
        cur_peak_min_pval <- signals$pval[i+1]
        cur_peak_top_snp <- signals$SNP[i+1]
        # and reset number of snps in peak
        cur_peak_n_snps <- 1
        i <- i + 1
      }
    }
    # there are other sig snps in LD with this SNP
    else
    {
      # update end of the current peak to the furthest snp in ld
      cur_peak_end <- cur_ld$snp2_bp[n_snps_in_ld]
      # update number of snps in peak
      cur_peak_n_snps <- cur_peak_n_snps + n_snps_in_ld
      # update cur_peak top snp info
        # get p-values of snps in ld with this snp
        # check if any have a p-val less than the current peak snp
        # if so, update the top snp info
      snps_in_ld <- filter(signals, SNP %in% c(cur_ld$snp2, signals$SNP[i]))
      if(min(snps_in_ld$pval, na.rm = TRUE) < cur_peak_min_pval)
      {
        cur_peak_min_pval <- min(snps_in_ld$pval, na.rm = TRUE)
        cur_peak_top_snp <- snps_in_ld$SNP[which.min(snps_in_ld$pval)]
      }
      # update i to skip to the furthest snp in ld with the current snp
      next_snp_name <- cur_ld$snp2[n_snps_in_ld]
      i <- which(signals$SNP==next_snp_name)
    }
  }
  if(nrow(peak_df) < 1)
  {
    print(str_c('No peaks found on chromosome ', chr, ' for ', trait))
    return()
  }
  return(peak_df)
}

# OPERATE ON SPECIFIED FILES
traits <- read_tsv(traitlist, col_names = 'trait')$trait
signals <- read_csv(signalfile) %>% 
	filter(CHROM==chrom)
ld <- read_csv(ldfile)

peaks <- tibble()

for(t in traits)
{
	cur_signals <- filter(signals, modelID==t)
	cur_peaks <- callLocalPeaks(cur_signals, ld, chrom, t)
	if(!is.null(cur_peaks))
	{
		peaks <- bind_rows(peaks, cur_peaks)
	}
}

write_csv(peaks, str_replace(ldfile, '.csv', str_c('_', tag, '_peaks.csv')), quote = 'needed')
