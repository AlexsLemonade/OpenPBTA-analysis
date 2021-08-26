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
  plot_title <- paste0("Telomerase activity scores for ", broad_histology)
  plot_fname <- file.path(plots_dir, paste0('EXTENDScores_', gsub(' ', '_',broad_histology), '.png'))
  table_fname <- file.path(results_dir, paste0('EXTENDScores_', gsub(' ', '_',broad_histology), '.tsv'))
  
  # create labels: count of samples per molecular subtype
  scores_mat <- scores_mat %>%
    mutate(molecular_subtype = paste0(molecular_subtype,' (N=',n_samples,')'))
  
  # create plot per molecular subtype
  if(nrow(scores_mat) > 1){
    
    # calculate adjusted p-value
    # default method = "kruskal.test" for multiple groups
    # adj <- compare_means(NormEXTENDScores ~ molecular_subtype, 
    #                      data = scores_mat, 
    #                      p.adjust.method = "bonferroni")
    
    # default method = "kruskal.test" for multiple groups
    # compare each group against "all"
    adj <- compare_means(NormEXTENDScores ~ molecular_subtype, 
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
    y_coord <-  max(scores_mat$NormEXTENDScores)
    
    # order of plot
    scores_mat$molecular_subtype <- factor(scores_mat$molecular_subtype, levels = scores_mat %>%
      group_by(molecular_subtype) %>%
      summarise(median = mean(NormEXTENDScores)) %>%
      arrange(median) %>%
      .$molecular_subtype)

    # plot
    p <- ggviolin(scores_mat, x = "molecular_subtype", y = "NormEXTENDScores", 
                   color = "molecular_subtype", 
                   palette = "simpsons", 
                   add = c("boxplot", "jitter"),
                   ggtheme = theme_pubr()) + 
      # manually add adjusted p-value from table generated above
      stat_pvalue_manual(adj, label = "p.signif", y.position = y_coord + 0.1)+
      # Add global p-value
      stat_compare_means(label.y = y_coord + 0.3) +
      xlab("Molecular subtype") +
      ylab("NormEXTENDScores") +
      rremove("legend") + 
      rotate_x_text(45)
    ggsave(plot = p, filename = plot_fname, width = 17, height = 8)
  }
}
