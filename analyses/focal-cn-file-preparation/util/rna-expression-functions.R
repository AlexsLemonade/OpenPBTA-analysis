# This script defines custom functions to be sourced in the 
# `rna-expression-validation.R` script of this module.
#
# Chante Bethell for CCDL 2019
#
# # #### USAGE
# This script is intended to be sourced in the 
# 'analyses/focal-cn-file-preparation/rna-expression-validation.R' script as 
# follows:
# 
# source(file.path("util", "rna-expression-functions.R"))

merge_expression <-
  function (copy_number_df, expression_df, metadata, filename){
    # Given the focal copy number data.frame already annotated with the
    # metadata, the RNA-seq expression data.frame, and the metadata, combine
    # the data into one data.frame and save data.frame as tsv file. 
    #
    # Args:
    #   copy_number_df: focal copy number data.frame annotated with the metadata
    #   expression_df: RNA-seq expression data.frame
    #   metadata: the relevant metadata data.frame
    #   filename: filename of the output tsv file
    #
    # Returns:
    #   combined_subset_df: data.frame with information from the focal CN, the
    #                       RNA-seq expression, and metadata data.frames, which
    #                       has been subsetted for plotting
    
    # Change expression data.frame to long tidy format
    long_expression <- expression_df %>%
      tidyr::gather(biospecimen_id, expression_value, -gene_id) %>%
      dplyr::distinct()
    long_expression$gene_id <-
      gsub(".*\\_", "", long_expression$gene_id)
    
    # Annotate expression data with metadata
    expression_metadata <- long_expression %>%
      dplyr::inner_join(metadata,
                        by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
      dplyr::select(
        biospecimen_id,
        expression_value,
        gene_id,
        Kids_First_Participant_ID,
        tumor_descriptor
      ) %>%
      dplyr::distinct()
    
    # Define mean and standard deviation for calculation of z-score
    mean <- mean(expression_metadata$expression_value)
    standard_deviation <- sd(expression_metadata$expression_value)
    
    # Merge Focal CN data.frame with RNA expression data.frame
    combined_df <- copy_number_df %>%
      dplyr::distinct() %>%
      dplyr::inner_join(
        expression_metadata,
        by = c(
          "Kids_First_Participant_ID",
          "gene_symbol" = "gene_id",
          "tumor_descriptor"
        )
      ) %>%
      dplyr::rename(biospecimen_id_cn = biospecimen_id.x, 
                    biospecimen_id_expression = biospecimen_id.y) %>%
      dplyr::mutate(z_score = log2(expression_value + 1) - mean / standard_deviation)
    
    # Define cutoffs for z_score
    z_minimum <- min(combined_df$z_score)
    z_median <- median(combined_df$z_score)
    z_mean <- mean(combined_df$z_score)
    z_third_quantile <- quantile(combined_df$z_score, prob = 0.75)
    z_maximun <- max(combined_df$z_score)
    
    combined_df <- combined_df %>%
      dplyr::mutate(expression_class = dplyr::case_when(
        z_score >= z_minimum & z_score < 0 ~ "Min-to-0",
        z_score == 0 ~ "0",
        z_score > 0 & z_score <= z_median ~ "0-to-Median",
        z_score > z_median & z_score <= z_third_quantile ~ "Median-to-3rd-Quantile",
        z_score > z_third_quantile & z_score <= z_mean ~ "3rd-Quantile-to-Mean",
        z_score > z_mean ~ "Mean-to-Max"))
    
    # Save results
    readr::write_tsv(combined_df, file.path(results_dir, filename))
    
  }

plot_stacked_expression <- function (cn_expression_df, filename_a, filename_b) {
  # Given a data.frame with annotated CN and RNA expression data, produce 
  # stacked barplots.
  #
  # Args: 
  #   cn_expression_df: data.frame with annotated CN and RNA expression data 
  #                     produced using `merge_expression` custom function
  #   filename_a: name to save the first produced plot as
  #   filename_b: name to save the second produced plot as
  #
  # Returns: 
  #   cn_expression_df: the given data.frame now filtered for the wanted gene
  #                     symbols
  
  # Filter for wanted gene symbols 
  cn_expression_df <- cn_expression_df %>%
    dplyr::filter(gene_symbol %in% goi_list$V1)
  
  # Plot and save 
  cn_expression_plot <- ggplot2::ggplot(cn_expression_df,
                                        ggplot2::aes(x = status, 
                                                     fill = expression_class)) +
    ggplot2::geom_bar(position = ggplot2::position_fill(reverse = TRUE)) +
    ggplot2::ylab("Proportion of called genes")
  
  ggplot2::ggsave(file.path(plots_dir, filename_a),
                  cn_expression_plot)
  
  # Plot per gene and save
  cn_expression_plot_per_gene <-
    cn_expression_plot + ggplot2::facet_wrap( ~ gene_symbol)
  
  ggplot2::ggsave(file.path(plots_dir, filename_b),
                  cn_expression_plot_per_gene)
  
  return(cn_expression_df)
  
} 

plot_mean_expression <- function (cn_df_loss, cn_df_non_loss, filename) {
  # Given a data.frame with mean expression values for loss calls and a 
  # data.frame with mean expression values for non-loss calls, produce and save
  # a scatterplot displaying this information across genes.
  #
  # Args: 
  #   cn_df_loss: data.frame containing just CN datat with the loss value in 
  #               `status` 
  #   cn_df_non_loss : data.frame containing just CN with non-loss values in 
  #                   `status`
  #   filename: name to save the output plot as
  #
  # Returns: 
  #   mean_combined_plot: a plot displaying the correlation between loss and 
  #                       non-loss calls 
  
  mean_loss <- cn_df_loss %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(loss_log_expression = mean(log2(expression_value + 1)))
  
  mean_non_loss <- cn_df_non_loss %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(non_loss_log_expression = mean(log2(expression_value + 1)))
  
  mean_combined <- mean_loss %>%
    dplyr::inner_join(mean_non_loss, by = "gene_symbol")
  
  mean_combined_plot <-
    ggplot2::ggplot(
      mean_combined,
      ggplot2::aes(x = non_loss_log_expression,
                   y = loss_log_expression,
                   col = gene_symbol)
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 1)
  ggplot2::ggsave(file.path(plots_dir, filename),
                  mean_combined_plot)
  
}
