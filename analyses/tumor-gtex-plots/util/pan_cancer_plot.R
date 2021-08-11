# plot/table of each cohort + cancer_group or cancer_group
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

pan_cancer_plot <- function(expr_mat_gene, hist_file, map_file,
                            analysis_type = c("cohort_cancer_group_level", "cancer_group_level"), 
                            plots_dir, results_dir, plot_width, plot_height, meta_file){
  
  # subset histology to minimal columns
  hist_file <- hist_file %>%
    dplyr::select(Kids_First_Biospecimen_ID, cohort, sample_type, cancer_group)
  
  # expression matrix to long format
  expr_mat_gene <- expr_mat_gene %>%
    tidyr::gather("Kids_First_Biospecimen_ID", "tpm", -c("gene"))
  
  # combine with histology file
  expr_mat_gene <- expr_mat_gene %>%
    inner_join(hist_file, by = "Kids_First_Biospecimen_ID")
  
  # get gene name and tumor cohorts
  gene_name <- unique(expr_mat_gene$gene)
  tumor_cohort <- paste0(unique(expr_mat_gene$cohort), collapse = ", ")
  
  # format x-axis labels and filter to n >= 5
  if(analysis_type == "cohort_cancer_group_level"){
    expr_mat_gene <- expr_mat_gene %>%
      group_by(cohort, cancer_group) %>%
      mutate(n_samples = n()) %>%
      filter(n_samples >= 5) %>%
      mutate(x_labels = paste0(cancer_group, ", ", cohort,  " (N = ", n_samples, ")"))
    
    # cohort name and title
    cohort_name <- tumor_cohort
    title <- paste(gene_name, "Gene Expression across cohorts", sep = "\n")
  } else if(analysis_type == "cancer_group_level") {
    expr_mat_gene <- expr_mat_gene %>%
      group_by(cancer_group) %>%
      mutate(n_samples = n()) %>%
      filter(n_samples >= 5) %>%
      mutate(x_labels = paste0(cancer_group, " (N = ", n_samples, ")"))
    
    # cohort name and title
    cohort_name <- "all_cohorts"
    title <- paste(gene_name, "Gene Expression across cancers", sep = "\n")
  }
  
  # reorder by median tpm
  fcts <- sort(unique(expr_mat_gene$x_labels))
  expr_mat_gene$x_labels <- factor(expr_mat_gene$x_labels, levels = fcts)
  
  # create unique title and filenames
  tumor_cohort_fname <- paste0(unique(expr_mat_gene$cohort), collapse = "_")
  fname <- paste(gene_name, tumor_cohort_fname, "pan_cancer", analysis_type, sep = "_")
  plot_fname <- paste0(fname, '.png')
  table_fname <- paste0('pan_cancer_plots_', analysis_type, '.tsv')
  
  # data-frame for metadata output 
  meta_df <- data.frame(Gene_symbol = gene_name, 
                        plot_type = "pan_cancer", 
                        Dataset = cohort_name,
                        Disease = NA,
                        analysis_type = analysis_type, 
                        plot_fname = plot_fname,
                        table_fname = table_fname)
  meta_file <- file.path(results_dir, 'metadata.tsv')
  if(!file.exists(meta_file)){
    write.table(x = meta_df, file = meta_file, sep = "\t", row.names = F, quote = F)
  } else {
    write.table(x = meta_df, file = meta_file, sep = "\t", row.names = F, col.names = F, quote = F, append = TRUE)
  }
  
  # boxplot
  cols <- c("Tumor" = "grey80")
  output_plot <- ggplot(expr_mat_gene, aes(x = x_labels, y = tpm, fill = sample_type)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    ylab("TPM") + xlab("") +
    theme_Publication2(base_size = 12) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(title) +
    scale_fill_manual(values = cols) + guides(fill = "none")
  ggsave(plot = output_plot, filename = file.path(plots_dir, plot_fname), device = "png", width = plot_width, height = plot_height)
  
  # output table of gene, median and sd
  output_table <- expr_mat_gene %>%
    group_by(gene, cohort, x_labels) %>%
    summarise(mean = mean(tpm),
              median = median(tpm),
              sd = sqrt(var(tpm))) %>%
    mutate(mean = round(mean, digits = 2),
           median = round(median, digits = 2),
           sd = round(sd, digits = 2))
  
  # replace cohort with "all_cohorts" for cancer_group_level
  if(analysis_type == "cancer_group_level"){
    output_table <- output_table %>%
      mutate(cohort = "all_cohorts")
  }
  
  # for now add dummy values for all other columns
  output_table <- output_table %>%
    dplyr::rename(Gene_symbol = gene,
                  Dataset = cohort) %>%
    mutate(Disease = gsub(" [(].*|[,].*", "", x_labels)) %>%
    inner_join(map_file, by = c("Gene_symbol" = "gene_symbol")) %>%
    dplyr::rename(Gene_Ensembl_ID = ensg_id) %>%
    dplyr::select(Gene_symbol, Gene_Ensembl_ID, Dataset, Disease, 
                  x_labels, mean, median, sd)
  table_fname <- file.path(results_dir, table_fname)
  if(!file.exists(table_fname)){
    write.table(x = output_table, file = table_fname, sep = "\t", row.names = F, quote = F)
  } else {
    write.table(x = output_table, file = table_fname, sep = "\t", row.names = F, col.names = F, quote = F, append = TRUE)
  }
}