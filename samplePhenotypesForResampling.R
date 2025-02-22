library(tidyverse)
# infile: path of csv file of phenotype data to sample 
# genotype: name of column with genotype information as a string
# trait: name of phenotype column to create samples of, as a string
# prop_keep: proportion of observations to retain as non-NA in each sample
# n_samples: number of samples to generate
# samples_per_file: number of samples to write per output file
samplePhenotypesForResampling <- function(infile, genotype, trait, prop_keep = 0.9, n_samples = 100, samples_per_file = 25)
{
  df <- read.csv(infile) %>% 
    select(all_of(c(genotype, trait)))
  n_obs <- nrow(df)
  obs_set_na <- as.integer((1 - prop_keep)*n_obs)
  n_files <- ceiling(n_samples/samples_per_file)
  
  samples <- tibble('{genotype}' := df[[genotype]])
  for(i in 1:n_samples)
  {
   samples <- samples %>% 
     mutate('{trait}_{i}':= replace(.data[[trait]], sample.int(1:n_obs, obs_set_na), NA))
  }
  
  samples <- select(samples, !any_of(trait))
  
  for(i in 1:n_files)
  {
    end_sample <- i*samples_per_file
    start_sample <- end_sample - samples_per_file + 1
    if(i==n_files){end_sample <- n_samples}
    
    samples_select <- str_c(trait, start_sample:end_sample, sep = '_')
    df_write <- select(samples, all_of(c(genotype, samples_select)))
    write.csv(df_write, str_c(str_remove(infile, '.csv'), '_', trait, '_', i, '.csv'), quote = FALSE, row.names = FALSE)
  }
}

