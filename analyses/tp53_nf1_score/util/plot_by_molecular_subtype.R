suppressPackageStartupMessages({
  library(ggpubr)
  library(ggplot2)
  library(tidyverse)
  library(ggsci)
})

plot_by_molecular_subtype <- function(scores_mat, plots_dir, results_dir){
  
  # plot title
  broad_histology <- unique(scores_mat$broad_histology)
  print(broad_histology)
  print(length(unique(scores_mat$molecular_subtype)))
  plot_title <- paste0("TP53 activity scores for ", broad_histology)
  plot_fname <- file.path(plots_dir, paste0('tp53_scores_vs_molecular_subtype_', gsub(' ', '_',broad_histology), '.png'))
  table_fname <- file.path(results_dir, paste0('tp53_scores_vs_molecular_subtype_', gsub(' ', '_',broad_histology), '.tsv'))
  
  # create labels: count of samples per molecular subtype
  scores_mat <- scores_mat %>%
    mutate(molecular_subtype = paste0(molecular_subtype,' (N=',n_samples,')'))
  
  # create plot per molecular subtype
  if(nrow(scores_mat) > 1){
    
    # default method = "kruskal.test" for multiple groups
    # compare each group against "all"
    adj <- compare_means(tp53_score ~ molecular_subtype, 
                         data = as.data.frame(scores_mat), 
                         p.adjust.method = "bonferroni", 
                         ref.group = ".all.") %>%
      rename('variable' = '.y.') %>%
      mutate(group1 = "all") %>%
      select(variable, group1, group2, p.format, p.adj, p.signif, method)
    
    # add significance codes instead of actual p-values
    adj$p.signif[adj$p.adj <= 0.0001] <- '****'
    adj$p.signif[adj$p.adj <= 0.001 & adj$p.adj > 0.0001] <- '***'
    adj$p.signif[adj$p.adj <= 0.01 & adj$p.adj > 0.001] <- '**'
    adj$p.signif[adj$p.adj <= 0.05 & adj$p.adj > 0.01] <- '*'
    adj$p.signif[adj$p.adj > 0.05] <- ''
    write.table(adj, file = table_fname, quote = F, sep = "\t", row.names = F)
    
    # get max y-axis for position of p-value labels
    y_coord <-  max(scores_mat$tp53_score)
    
    # order of plot
    scores_mat$molecular_subtype <- factor(scores_mat$molecular_subtype, levels = scores_mat %>%
                                             group_by(molecular_subtype) %>%
                                             summarise(median = mean(tp53_score)) %>%
                                             arrange(median) %>%
                                             .$molecular_subtype)
    
    # set width for mol subtypes > 10
    if(length(unique(scores_mat$molecular_subtype)) > 10){
      width = 38
    } else {
      width = 17
    }
    
    # plot
    p <- ggviolin(scores_mat, x = "molecular_subtype", y = "tp53_score", 
                  color = "molecular_subtype", 
                  palette = "simpsons", 
                  add = c("boxplot", "jitter"),
                  ggtheme = theme_pubr()) + 
      # manually add adjusted p-value from table generated above
      stat_pvalue_manual(adj, label = "p.signif", y.position = y_coord + 0.3)+
      # Add global p-value
      stat_compare_means(label.y = y_coord + 0.5) +
      xlab("Molecular subtype") +
      ylab("TP53 scores") +
      rremove("legend") + 
      rotate_x_text(45)
    ggsave(plot = p, filename = plot_fname, width = width, height = 8)
  }
}