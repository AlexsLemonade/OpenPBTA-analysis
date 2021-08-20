suppressPackageStartupMessages({
  library(ggpubr)
  library(ggplot2)
  library(tidyverse)
})

theme_set(theme_classic() +
            theme(axis.text.x = element_text(angle = 50, size = 7, vjust = 1, hjust = 1),
                  legend.position = "top",
                  legend.key.size = unit(0.3, "cm"),
                  legend.key.width = unit(0.3, "cm"),
                  legend.title = element_text(size = 7),
                  legend.text = element_text(size = 6)
            )
)

boxplot_by_molecular_subtype <- function(scores_mat, output_dir){
  
  # plot title
  short_histology <- unique(scores_mat$short_histology)
  plot_title <- paste0("Telomerase activity scores for ", short_histology)
  fname <- file.path(output_dir, paste0('EXTENDScores_', gsub(' ', '_',short_histology), '.png'))
  
  # create labels: count of samples per molecular subtype
  scores_mat <- scores_mat %>%
    group_by(short_histology, molecular_subtype, NormEXTENDScores) %>%
    unique() %>%
    mutate(label = n()) %>%
    mutate(label = paste0(molecular_subtype,' (',label,')'))
  
  # create boxplot per molecular subtype
  if(nrow(scores_mat) > 1){
    y_coord <-  max(scores_mat$NormEXTENDScores) + 0.1
    p <- ggplot(scores_mat, aes(x = fct_reorder(molecular_subtype, NormEXTENDScores, .desc = TRUE),
                           y = NormEXTENDScores)) + 
      geom_boxplot(width = 0.5, size = 0.1, notch = FALSE, outlier.size = 0, outlier.shape = NA, fill = "pink") + 
      geom_jitter(shape = 16, width = 0.1, size = 0.4) +
      # pairwise comparison against all
      stat_compare_means(label = "p.signif", method = "t.test", label.y = y_coord, ref.group = ".all.", color = "darkred") + 
      stat_compare_means(method = "anova", label.y = y_coord + 0.1, color = "darkred", label.y.npc = "top", label.x.npc = "middle") + 
      xlab('') + ggtitle(plot_title)
    ggsave(plot = p, filename = fname, width = 7, height = 5)
  }
}
