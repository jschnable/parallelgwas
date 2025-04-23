library(tidyverse)

summarisePeaks <- function(path)
{
  files <- Sys.glob(path)
  peaks <- tibble()
  for(f in files)
  {
    df <- read_csv(f) %>% 
      mutate(filename = f)
    peaks <- bind_rows(peaks, df)
  }
  return(peaks)
}
