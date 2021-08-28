# function to plot ROC curve
# Author: Jo Lynne Rokita

suppressPackageStartupMessages({
  library(ggpubr)
  library(ggplot2)
})

plot_roc <- function(roc_file, plots_dir, fname){
  roc_df <- readr::read_tsv(roc_file,
                            col_types = readr::cols(
                              .default = readr::col_double(),
                              gene = readr::col_character(),
                              shuffled = readr::col_character()
                            )
  )
  
  curve_colors <- c(
    "#313695",
    "#35978f",
    "#c51b7d",
    "#313695",
    "#35978f",
    "#c51b7d"
  )
  
  curve_labels <- c(
    "TP53 False" = "TP53 (AUROC = 0.89)",
    "Ras False" = "Ras (AUROC = 0.55)",
    "NF1 False" = "NF1 (AUROC = 0.77)",
    "TP53 True" = "TP53 Shuffled (AUROC = 0.46)",
    "Ras True" = "Ras Shuffled (AUROC = 0.55)",
    "NF1 True" = "NF1 Shuffled (AUROC = 0.43)"
  )
  roc_df$model_groups <- paste(roc_df$gene, roc_df$shuffled)
  roc_df$model_groups <- factor(roc_df$model_groups, levels = names(curve_labels))
  
  p <- ggplot(roc_df,
              aes(x = fpr,
                  y = tpr,
                  color = model_groups)) +
    coord_fixed() +
    geom_step(aes(linetype = shuffled), size = 0.7) +
    geom_segment(aes(x = 0,
                     y = 0,
                     xend = 1,
                     yend = 1),
                 linetype = "dashed",
                 color = "black") +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    scale_color_manual(name = "Classifier",
                       values = curve_colors,
                       labels = curve_labels) +
    scale_linetype_manual(name = "Data",
                          values = c("solid",
                                     "dashed"),
                          labels = c("True" = "Shuffled",
                                     "False" = "Real")) +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") + theme_pubr()
  ggsave(p, filename = file.path(plots_dir, fname))
}