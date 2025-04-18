library(tidyverse)
library(scales)
# data: dataframe containing chromosome, bp, and significance (either RMIP or p-values) values for each marker
# sig: column containing the significance values
# multitrait: logical indicating if the plot should represent GWAS results for multiple phenotypes
# trait: column indicating which trait the marker significance value is for if multitrait is TRUE
# resampling: logical indicating if the data is from a resampling-based GWAS or not
# chromLengths: vector of number of basepairs in each chromosome; if NULL, this is estimated based on markers in data
# chr: name of column with chromosome numbers 
# bp: name of column with basepair number for each marker within its chromosome
# threshold: position of significance threshold to plot; default is 0.1, assuming the default of resampling = TRUE
# main: string to plot as plot main title
# colors: vector of hex codes or strings of color names in the R color list to use; if provided, must be the length of either the number of chromosomes (if multitrait=FALSE) or the number of traits (if multitrait = TRUE)
# theme: object containing theme info to use instead of theme elements hard-coded here
# species: name of species; currently supported: 'maize', 'sorghum'
plotManhattan <- function(data, sig, multitrait=FALSE, trait=NULL, resampling=TRUE, chromLengths=NULL, chr=CHROM, bp=POS, threshold=0.1, main=NULL, colors=NULL, theme=NULL, species)
{
  ylab <- 'RMIP'
  theme_use <- theme

  if(species=='maize' & is.null(chromLengths)){chromLength <-  tibble(max_bp = c(308452471, 243675191, 238017767, 250330460, 226353449, 181357234, 185808916, 182411202, 163004744, 152435371), {{ chr }} := 1:10)}
  if(species=='sorghum' & is.null(chromLengths)){chromLength <- tibble(max_bp = c(85112863, 79114963, 80873341, 71215609, 77058072, 62713908, 68911884, 65779274, 63277606, 62870657), {{ chr }} := 1:10)}
  
  if(is.null(theme))
  {
    theme_use <- theme_minimal() +
      theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                       vjust = 0.5, hjust = 0.5),
            axis.text.y = element_text(size = 11, color = 'black', vjust = 0, hjust = 0.5),
            legend.text = element_text(size = 11, color = 'black'),
            plot.title = element_text(size = 11, color = 'black', vjust = 0, hjust = 0.5),
            text = element_text(size = 11, color = 'black'),
            legend.position = 'none',
            line = element_line(color = 'black', linewidth = 1),
            panel.grid = element_blank())
  }
  
  data_cum <- data %>%
    group_by( {{ chr }}) %>%
    summarise(max_bp = max({{ bp }}, na.rm = TRUE)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(c({{ chr }}, bp_add))
  
  if(!is.null(chromLengths))
  {
    data_cum <- data %>%
      group_by({{ chr }}) %>%
      summarise(n_snps = n()) %>%
      full_join(chromLengths, join_by('{chr}')) %>%
      arrange({{ chr }}) %>% 
      mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
      select(c({{ chr }}, bp_add))
  }
  
  df <- data %>% 
    inner_join(data_cum, join_by({{ chr }})) %>%
    rowwise() %>%
    mutate(loc = {{ bp }} + bp_add)

   if(!resampling)
   {
     ylab <-  expression(-log[10]~'(p)')
     df <- df %>% 
       rowwise() %>%
       mutate({{ sig }} := -1*log10({{ sig }}))
   }
  
  x_axis_set <- df %>% 
    group_by({{ chr }}) %>% 
    summarise(center = mean(loc, na.rm = TRUE))
  
  if(multitrait)
  {
    if(is.null(colors))
    {
      n_traits <- df %>% 
        group_by({{ trait}}) %>%
        summarise(n = n())
      colors <- hue_pal()(dim(n_traits)[1])
    }
    manhattan <- ggplot(df, aes(loc, {{ sig }}, color = as.factor(.data[[deparse(substitute(trait))]]))) + 
      geom_point() + 
      geom_hline(yintercept = threshold, linetype = 2, color = 'black') +
      scale_x_continuous(label = x_axis_set[[deparse(substitute(chr))]], 
                         breaks = x_axis_set$center) +
      scale_color_manual(values = colors) + 
      labs(title = main, x = 'Chromosome', y = ylab, color = NULL) + 
      theme_use + 
      theme(legend.position = 'top')
    print(manhattan)
    return(manhattan)
  }
  else
  {
    if(is.null(colors))
    {
      colors <- hue_pal()(dim(data_cum)[1])
    }
    manhattan <- ggplot(df, aes(loc, {{ sig }}, color = as.factor({{ chr }}))) + 
      geom_point() + 
      geom_hline(yintercept = threshold, linetype = 2, color = 'black') +
      scale_x_continuous(label = x_axis_set[[deparse(substitute(chr))]], 
                         breaks = x_axis_set$center) +
      scale_color_manual(values = colors) + 
      labs(title = main, x = 'Chromosome', y = ylab, color = NULL) + 
      theme_use
    print(manhattan)
    return(manhattan)
  }
}
