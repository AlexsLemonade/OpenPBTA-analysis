plot_by_cancer_predisposition <- function(lfs, plot_dir){
  cancer_predisposition <- unique(lfs$cancer_predispositions)
  fname <- paste0('tp53_scores_vs_tp53_altered_status_', cancer_predisposition, '.png')
  fname <- file.path(plot_dir, fname)
  
  # use kruskal test for > 2 levels
  if(length(unique(lfs$tp53_altered)) > 2){
    method <- "kruskal.test"
  } else {
    method <- "wilcox.test"
  }
  p <- ggviolin(lfs, x = "tp53_altered", y = "tp53_score", 
           color = "tp53_altered", 
           palette = "jco",
           order = c("activated", "loss", "other"),
           add = c("boxplot", "jitter"),  
           ggtheme = theme_pubr()) +
    # Add pairwise comparisons p-value
    stat_compare_means(method = method, label.y = 1.2, label.x.npc = "center") +
    xlab("TP53 altered status") +
    ylab("TP53 score") +
    rremove("legend")
  ggsave(plot = p, filename = fname)
}
