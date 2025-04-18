library(tidyverse)
# data: Data frame to split into separate files by col
# cols: cols to split into separate files; tidy-select
# out: path to write file subsets to as csv files
# subsetSize: number of columns of cols per out file
splitDataFrame <- function(data, cols, out, subsetSize)
{
  metadata <- select(data, !cols)
  data_split <- select(data, cols)
  cols_to_split <- ncol(data_split)
  n_subsets <- ceiling(cols_to_split/subsetSize)
  
  for(i in 1:n_subsets)
  {
    outfile <- paste0(out, i, '.csv')
    end_index <- i*subsetSize
    start_index <- end_index - (subsetSize - 1)
    if(end_index > cols_to_split){end_index <- cols_to_split}
    subset <- data_split[, start_index:end_index]
    subset <- bind_cols(metadata, subset)
    write.csv(subset, outfile, quote = FALSE, row.names = FALSE)
  }
}
