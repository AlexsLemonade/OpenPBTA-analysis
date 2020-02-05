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

calculate_z_score <- function (expression_data, rnaseq_ind) {
  # Given an expression data matrix, filter for independent samples and
  # calculate the log2 expression values and z-scores.
  #
  # Args:
  #   expression_data: matrix with RNA expression data
  #
  # Return:
  #   long_expression: expression data.frame filtered for independent samples
  #                    and now in long format, and with the added columns:
  #                    `z_score`, `expression_class`
  
  # log2 transform expression matrix
  log2_exp <- log2(expression_data + 1)
  
  # z-score log2 transformed expression matrix
  zscore_expression <- scale(t(log2_exp), 
                             center = TRUE,
                             scale = TRUE)
  
  # Change z-scored expression data.frame to long tidy format
  long_expression <- zscore_expression %>%
    reshape2::melt() %>%
    dplyr::rename(
      biospecimen_id = Var1,
      gene_id = Var2,
      z_score = value
    ) %>%
    dplyr::filter(biospecimen_id %in% rnaseq_ind) %>%
    dplyr::mutate(
      gene_id = gsub(".*\\_", "", gene_id),
      # Trim the gene ids to
      # include only the gene symbol
      expression_class = dplyr::case_when(
        z_score < -2 ~ "z < -2",
        z_score < -1 ~ "-2 ≤ z < -1",
        z_score < -0.5 ~ "-1 ≤ z < -0.5",
        z_score < 0 ~ "-0.5 ≤ z < 0",
        z_score < 0.5 ~ "0 ≤ z < 0.5",
        z_score < 1 ~ "0.5 ≤ z < 1",
        z_score < 2 ~ "1 ≤ z < 2",!(is.na(z_score)) ~ "z ≥ 2"
      ) %>%
        ordered(
          levels = c(
            "z < -2",
            "-2 ≤ z < -1",
            "-1 ≤ z < -0.5",
            "-0.5 ≤ z < 0",
            "0 ≤ z < 0.5",
            "0.5 ≤ z < 1",
            "1 ≤ z < 2",
            "z ≥ 2"
          )
        )
    ) %>%
    dplyr::distinct()
  
}

merge_expression <-
  function (copy_number_df,
            expression_df,
            metadata,
            filename_lead) {
    # Given the focal copy number data.frame already annotated with the
    # metadata, the RNA-seq expression data.frame, and the metadata, combine
    # the data into one data.frame and save data.frame as tsv file.
    #
    # Args:
    #   copy_number_df: focal copy number data.frame
    #   expression_df: RNA-seq expression data.frame
    #   metadata: the relevant metadata data.frame
    #   filename_lead: the lead for the filename of the output tsv file
    #
    # Returns:
    #   combined_df: data.frame with information from the focal CN, the
    #                       RNA-seq expression, and metadata data.frames
    
    # Annotate expression data with metadata
    expression_metadata <- expression_df %>%
      dplyr::inner_join(metadata,
                        by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
      dplyr::select(
        gene_id,
        Kids_First_Participant_ID,
        sample_id,
        biospecimen_id,
        z_score,
        expression_class,
        tumor_descriptor
      ) %>%
      dplyr::distinct()
    
    # Annotate focal CN data with metadata
    cn_metadata <- copy_number_df %>%
      dplyr::inner_join(metadata,
                        by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
      dplyr::select(
        gene_symbol,
        Kids_First_Participant_ID,
        sample_id,
        biospecimen_id,
        status,
        copy_number,
        tumor_ploidy,
        tumor_descriptor
      ) %>%
      dplyr::distinct()
    
    # Merge Focal CN data.frame with RNA expression data.frame
    combined_df <- expression_metadata %>%
      dplyr::left_join(
        cn_metadata,
        by = c("sample_id",
               "gene_id" = "gene_symbol",
               "tumor_descriptor"),
        suffix = c("_expression", "_cn")
      )
    
    # Save results
    readr::write_tsv(combined_df, file.path(results_dir, paste0(filename_lead, "_expression_merged.tsv")))
    
    return(combined_df)
    
  }

plot_stacked_expression <- function (cn_expression_loss_df,
                                     cn_expression_neutral_df,
                                     cn_expression_zero_df,
                                     plotname_lead) {
  # Given a data.frame with annotated CN and RNA expression data, produce a
  # stacked barplot for loss calls, neutral calls, and instances where
  # `copy_number` == 0.
  #
  # Args:
  #   cn_expression_loss_df: data.frame with annotated CN and RNA expression
  #                          data produced using `merge_expression` custom
  #                          function and filtered for loss calls
  #   cn_expression_neutral_df: data.frame with annotated CN and RNA expression
  #                             data produced using `merge_expression` custom
  #                             function and filtered for neutral calls
  #   cn_expression_zero_df: data.frame with annotated CN and RNA expression
  #                          data produced using `merge_expression` custom
  #                          function and filtered for `copy_number` = 0
  #   plotname_lead: lead for name to save the combined stacked barplot as
  #
  # Returns:
  #   Saves stacked barplot.
  
  # Bind input data.frame rows
  cn_expression_df <- dplyr::bind_rows(
    loss = cn_expression_loss_df,
    zero = cn_expression_zero_df,
    neutral = cn_expression_neutral_df,
    .id = "status_name"
  )

  # Plot and save
  cn_expression_plot <- ggplot2::ggplot(cn_expression_df,
                                        ggplot2::aes(x = status_name,
                                                     fill = expression_class)) +
    ggplot2::geom_bar(position = ggplot2::position_fill(reverse = TRUE)) +
    ggplot2::ylab("Proportion of called genes") +
    ggplot2::labs(title = toupper(paste0(plotname_lead, "_stacked_plot")))

  ggplot2::ggsave(file.path(plots_dir, paste0(plotname_lead, "_stacked_plot.png")),
                  cn_expression_plot)

}

plot_mean_expression <- function (cn_expression_loss_df,
                                  cn_expression_neutral_df,
                                  cn_expression_zero_df,
                                  plotname_lead) {
  # Given a data.frame with expression values for all CN calls and a
  # data.frame with expression values for loss calls, produce and save
  # a scatterplot displaying the correlation of loss and neutral calls across
  # genes.
  #
  # Args:
  #   cn_expression_loss_df: data.frame with annotated CN and RNA expression
  #                          data produced using `merge_expression` custom
  #                          function and filtered for loss calls
  #   cn_expression_neutral_df: data.frame with annotated CN and RNA expression
  #                             data produced using `merge_expression` custom
  #                             function and filtered for neutral calls
  #   cn_expression_zero_df: data.frame with annotated CN and RNA expression
  #                          data produced using `merge_expression` custom
  #                          function and filtered for `copy_number` = 0
  #   plotname_lead: name to save the output neutral/loss correlation plot
  #
  # Returns:
  #   The above named plots are saved in `plots_dir`
  
  # Calculate the mean of log2 expression values
  mean_loss <- cn_expression_loss_df %>%
    dplyr::group_by(gene_id, is_driver_gene) %>%
    dplyr::summarise(mean_loss_zscore = mean(z_score))
  
  mean_neutral <- cn_expression_neutral_df %>%
    dplyr::group_by(gene_id, is_driver_gene) %>%
    dplyr::summarise(mean_neutral_zscore = mean(z_score))
  
  mean_zero <- cn_expression_zero_df %>%
    dplyr::group_by(gene_id, is_driver_gene) %>%
    dplyr::summarise(mean_zero_zscore = mean(z_score))
  
  # Combine the mean values to be plotted
  mean_combined_loss_neutral <- mean_loss %>%
    dplyr::inner_join(mean_neutral, by = c("gene_id", "is_driver_gene"))
  
  mean_combined_zero_neutral <- mean_zero %>%
    dplyr::inner_join(mean_neutral, by = c("gene_id", "is_driver_gene"))
  
  # Plot neutral/loss mean values
  mean_combined_plot_loss <-
    ggplot2::ggplot(
      mean_combined_loss_neutral,
      ggplot2::aes(x = mean_neutral_zscore,
                   y = mean_loss_zscore,
                   col = is_driver_gene)
    ) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::geom_abline() +
    ggplot2::labs(title = toupper(paste0(plotname_lead, "_loss_cor_plot")))
  ggplot2::ggsave(file.path(plots_dir, paste0(plotname_lead, "_loss_cor_plot.png")),
                  mean_combined_plot_loss)
  
  # Plot neutral/zero mean values
  mean_combined_plot_zero <-
    ggplot2::ggplot(
      mean_combined_zero_neutral,
      ggplot2::aes(x = mean_neutral_zscore,
                   y = mean_zero_zscore,
                   col = is_driver_gene)
    ) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::geom_abline() +
    ggplot2::labs(title = toupper(paste0(plotname_lead, "_zero_cor_plot")))
  ggplot2::ggsave(file.path(plots_dir, paste0(plotname_lead, "_zero_cor_plot.png")),
                  mean_combined_plot_zero)
  
}
