# base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git")) 
analysis_dir <- file.path(root_dir, "analyses", "tp53_nf1_score")

# source function to plot ROC curve
source(file.path(analysis_dir, "util", "plot_roc.R"))

# output directories
plots_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results")

# read data obtained from 06-evaluate-classifier.py
roc_file_polya <- readr::read_tsv(file.path(results_dir, "polya_TP53_roc_threshold_results.tsv"))
roc_file_polya_shuff <- readr::read_tsv(file.path(results_dir, "polya_TP53_roc_threshold_results_shuffled.tsv"))

roc_file_stranded <- readr::read_tsv(file.path(results_dir, "stranded_TP53_roc_threshold_results.tsv"))
roc_file_stranded_shuff <- readr::read_tsv(file.path(results_dir, "stranded_TP53_roc_threshold_results_shuffled.tsv"))

# call function to plot ROC curve
plot_roc(roc_df =  roc_file_polya %>%
           rbind(roc_file_polya_shuff), plots_dir, fname = "polya_TP53_roc.png")
plot_roc(roc_df =  roc_file_stranded %>%
           rbind(roc_file_stranded_shuff), plots_dir, fname = "stranded_TP53_roc.png")

