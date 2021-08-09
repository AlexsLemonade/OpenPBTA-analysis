# plot/table of each cohort + cancer_group vs GTEx subgroups
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

tumor_normal_gtex_plot <- function(expr_mat_gene, hist_file, map_file,
                                   analysis_type = c("cohort_cancer_group_level", "cancer_group_level"), 
                                   plots_dir, results_dir, plot_width, plot_height, meta_file){
  
  # standardize groups:
  # create a single group variable for both cancer_group and gtex_subgroup 
  hist_file <- hist_file %>%
    mutate(group = ifelse(sample_type == "Normal", gtex_subgroup, cancer_group)) %>%
    dplyr::select(Kids_First_Biospecimen_ID, cohort, sample_type, gtex_subgroup, cancer_group, group)
  
  # expression matrix to long format
  expr_mat_gene <- expr_mat_gene %>%
    tidyr::gather("Kids_First_Biospecimen_ID", "tpm", -c("gene"))
  
  # combine with histology file
  expr_mat_gene <- expr_mat_gene %>%
    inner_join(hist_file, by = "Kids_First_Biospecimen_ID")
  
  # format x-axis labels and filter to n >= 5
  if(analysis_type == "cohort_cancer_group_level"){
    expr_mat_gene <- expr_mat_gene %>%
      group_by(cohort, group) %>%
      mutate(n_samples = n()) %>%
      filter(n_samples >= 5) %>%
      mutate(x_labels = paste0(group, " (N = ", n_samples, ")"))
  } else {
    expr_mat_gene <- expr_mat_gene %>%
      group_by(group) %>%
      mutate(n_samples = n()) %>%
      filter(n_samples >= 5) %>%
      mutate(x_labels = paste0(group, " (N = ", n_samples, ")"))
  }
  
  # get a vector of all tumor groups
  cohort_cancer_groups <- expr_mat_gene %>%
    filter(sample_type == "Tumor") %>%
    .$x_labels %>%
    unique()
  
  # loop through cancer + cohort groups
  for(i in 1:length(cohort_cancer_groups)){
    
    # subset to specific cohort_cancer_group and all gtex  
    expr_mat_gene_subset <- expr_mat_gene %>%
      filter(x_labels %in% cohort_cancer_groups[i] | cohort == "GTEx")
    
    # reorder by alphabet 
    fcts <- expr_mat_gene_subset %>% 
      dplyr::select(sample_type, group, x_labels) %>%
      unique() %>%
      arrange(desc(sample_type), group) %>%
      .$x_labels 
    expr_mat_gene_subset$x_labels <- factor(expr_mat_gene_subset$x_labels, levels = fcts)
    
    # create unique title and filenames
    gene_name <- unique(expr_mat_gene_subset$gene)
    # cohorts <- paste0(unique(expr_mat_gene_subset$cohort), collapse = ", ")
    tumor_cohort <- expr_mat_gene_subset %>%
      filter(sample_type == "Tumor") %>%
      .$cohort %>% 
      unique() %>%
      paste0(collapse = ", ")
    if(analysis_type == "cohort_cancer_group_level"){
      cohort_name <- tumor_cohort
    } else {
      cohort_name <- "all_cohorts"
    }
    cancer_group_name <- expr_mat_gene_subset %>%
      filter(sample_type == "Tumor") %>%
      .$cancer_group %>%
      unique()
    # replace semi-colon, forward slash and spaces with hyphen in filenames
    cancer_group_name_fname <- gsub('/| |;', '-', cancer_group_name) 
    
    # create title and filename prefix
    if(analysis_type == "cohort_cancer_group_level"){
      title <- paste(gene_name, 
                     paste(tumor_cohort, cancer_group_name, "vs. GTEx", sep = " "), sep = "\n")
      fname <- paste(gene_name, tumor_cohort, cancer_group_name_fname, "vs_GTEx", analysis_type, sep = "_")
    } else {
      title <- paste(gene_name,
                     paste(cancer_group_name, "vs. GTEx", sep = " "), sep = "\n")
      fname <- paste(gene_name, cancer_group_name_fname, "vs_GTEx", analysis_type, sep = "_")
    }
    plot_fname <- paste0(fname, '.png')
    table_fname <- paste0('tumor_normal_gtex_plots_', analysis_type, '.tsv')
    
    # data-frame for metadata output 
    meta_df <- data.frame(Gene_symbol = gene_name, 
                          plot_type = "tumor_normal_gtex", 
                          Dataset = cohort_name,
                          Disease = gsub(" [(].*|[,].*", "", cohort_cancer_groups[i]),
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
    cols <- c("Normal" = "grey80", "Tumor" = "red3")
    output_plot <- ggplot(expr_mat_gene_subset, aes(x = x_labels, y = tpm, fill = sample_type)) +
      stat_boxplot(geom ='errorbar', width = 0.2) +
      geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
      ylab("TPM") + xlab("") +
      theme_Publication2(base_size = 12) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggtitle(title) +
      scale_fill_manual(values = cols) + guides(fill = "none")
    ggsave(plot = output_plot, filename = file.path(plots_dir, plot_fname), device = "png", width = plot_width, height = plot_height)
    
    # output table of gene, median and sd
    output_table <- expr_mat_gene_subset %>%
      group_by(gene, x_labels) %>%
      summarise(mean = mean(tpm),
                median = median(tpm),
                sd = sqrt(var(tpm))) %>%
      mutate(mean = round(mean, digits = 2),
             median = round(median, digits = 2),
             sd = round(sd, digits = 2))
    
    # for now add dummy values for all other columns
    output_table <- output_table %>%
      dplyr::rename(Gene_symbol = gene) %>%
      mutate(Dataset = cohort_name, 
             Disease = gsub(" [(].*|[,].*", "", cohort_cancer_groups[i])) %>%
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
}