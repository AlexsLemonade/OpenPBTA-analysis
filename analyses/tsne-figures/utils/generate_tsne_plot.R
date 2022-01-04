
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

generate_tsne_plot <- function(df, cancer_group, col_list){
  # get all the subtypes 
  subtypes <- df %>% 
    filter(subtypes_to_plot != "Other CNS tumor") %>% 
    pull(subtypes_to_plot) %>% 
    unique()
  # reorder the subtypes
  df$subtypes_to_plot <- factor(df$subtypes_to_plot, levels=c(subtypes, "Other CNS tumor"))
  
  # get the plots
  plot <- df %>%
    ggplot(aes(x = X1, 
               y = X2,
               color = subtypes_to_plot)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_bw() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_color_manual(values=c(sample(col_list, length(subtypes)),
                                "#bebebe")) 
  
  # write out the plots
  pdf(file.path(plots_dir, paste0("umap-", cancer_group, "-subtypes.pdf")), width = 11, height = 8.5)
  print(plot)
  dev.off()
  
  # return the plot
  return(plot)
}
