# function to plot ROC curve
# Author: Jo Lynne Rokita

suppressPackageStartupMessages({
  library(ggpubr)
  library(ggplot2)
  library(ggsci)
})

plot_roc <- function(roc_df, plots_dir, fname){

  # add legend 
  roc_df <- roc_df %>%
    mutate(auroc = round(auroc, 2),
           shuffled = ifelse(shuffled, 'TP53 Shuffle', 'TP53'),
           Classifier = paste0(shuffled, ' (AUROC = ', auroc,')'))

  # plot 
  p <- ggplot(roc_df, aes(x = fpr,  y = tpr)) +
    coord_fixed() +
    geom_step(aes(linetype = Classifier, color = Classifier), size = 0.7) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "black") +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") + theme_pubr() + scale_color_simpsons()
  ggsave(p, filename = file.path(plots_dir, fname), width = 8, height = 8)
}
