# Functions for making mutational signature plots and calculations
#
# C. Savonen for ALSF - CCDL
# 2019
#
################################################################################
sample_mut_sig_plot <- function(which_sig_list, label = "none", output_dir = getwd(),
                                sample_ids) {
  # Given a list of `deconstructSigs::WhichSignature` output, make plots for each sample using
  # deconstructSigs::plotSignatures
  #
  # Args:
  #
  #   which_sig_list: a list of `whichSignature` output for each sample
  #   label: the label you would like associated with the plot as a character string.
  #   sample_ids: vector of the sample ids.
  #   output_dir: where the plots should be saved. If this directory doesn' exist,
  #               will create it. Default is current directory.
  # Returns:
  #   Saved png mutation signature plots to the specified directory.

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  for (sample_num in 1:length(which_sig_list)) {
    png(file.path(output_dir, paste0(sample_ids[sample_num], "_", label, "_mutation_sig.png")))
    plotSignatures(which_sig_list[[sample_num]], sub = sample_ids[sample_num])
    dev.off()
  }
}

calc_mut_per_sig <- function(which_sig_list,
                             muts_per_sample,
                             wgs_genome_size,
                             wxs_exome_size,
                             maf_df) {
  # Given a list of `deconstructSigs::whichSignature` output, calculate the
  # mutations per signature per Mb of the genome/exome.
  #
  # Args:
  #
  #   which_sig_list: a list of `whichSignature` output
  #   muts_per_sample: a vector with the total mutation signature counts for each sample
  #   wgs_genome_size: size of the WGS genome in bp
  #   wxs_genome_size: size of the WXS exome in bp
  #   maf_df: a data.frame with `short_histology` and `experimental strategy` information columns
  #
  # Returns:
  #   A data.frame that is samples x signatures and has the number of mutations
  #   per mb for each signature x sample combo

  # Count the total number of signature mutations for each sample
  total_muts <- apply(sigs_input, 1, sum)

  # Pull out the signature weights and make into matrix
  sig_weights <- do.call(
    "rbind.data.frame",
    lapply(which_sig_list, function(sample_data) sample_data$weights)
  ) %>%
    tibble::rownames_to_column("Tumor_Sample_Barcode")

  # Calculate the number of mutations contributing to each signature
  # Here the weight is multiplied by the total number of signature mutations. 
  sig_num_df <- apply(sig_weights[, -1], 2, function(row_weights) row_weights * total_muts) %>%
    as.data.frame() %>%

    # Add back the sample ids
    dplyr::mutate(Tumor_Sample_Barcode = sig_weights$Tumor_Sample_Barcode) %>%

    # Join the short_histology and experimental stategy information
    dplyr::inner_join(dplyr::select(
      maf_df,
      "Tumor_Sample_Barcode",
      "short_histology",
      "experimental_strategy"
    ),
    by = "Tumor_Sample_Barcode"
    ) %>%
    # Get rid of Panel samples
    dplyr::filter(experimental_strategy != "Panel") %>%

    # Reformat for plotting
    reshape2::melt() %>%

    # Only keep distinct
    dplyr::distinct(Tumor_Sample_Barcode, variable, .keep_all = TRUE) %>%

    # Add genome size and calculate the mutation per this column
    dplyr::mutate(
      genome_size = dplyr::recode(experimental_strategy,
        "WGS" = wgs_size,
        "WXS" = wxs_size
      ),
      mut_per_mb = value / (genome_size / 10^6)
    )
  
  return(sig_num_df)
}

bubble_matrix_plot <- function(sig_num_df, label = "none") {
  # Given the data.frame output from `calc_mut_per_sig` that has the number of 
  # mutations per Mb that belong to each sample x signature combo, make a bubble 
  # matrix plot that displays this information for all samples but summarized by 
  # histology.
  #
  # Args:
  #
  #   sig_num_df: a data.frame with number of mutations per Mb that belong to
  #               each sample x signature combo
  #   label: a character string for the title of the plot to be passed to ggplot2::ggtitle
  #
  # Returns:
  #   A bubble matrix plot with the number of mutations per Mb and proportion of 
  #   tumors with a non-zero weight for all samples summarized by histology.
  #
  # Summarize the mut_per_mb column by histology
  grouped_sig_num <- sig_num_df %>%
    dplyr::arrange(variable) %>%
    dplyr::group_by(short_histology, variable) %>%
    dplyr::summarize(
      # Calculate the proportion of tumors with a non-zero weight for each signature
      prop_tumors = sum(value > 0) / length(value),
      # Calculate the median number of mutations per Md
      med_num = median(mut_per_mb)
    )

  # Make the bubble matrix plot
  ggplot2::ggplot(grouped_sig_num, ggplot2::aes(x = short_histology, y = forcats::fct_rev(variable))) +
    ggplot2::geom_point(ggplot2::aes(color = med_num, size = prop_tumors)) +
    ggplot2::scale_size("Proportion of Samples", range = c(0, 4)) +
    ggplot2::scale_color_distiller("Median Number of Mutations per Mb", palette = "YlGnBu") +
    ggplot2::theme_classic() +
    # Make labels on top
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 55, hjust = -.01),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(.75)),
      legend.key.size = ggplot2::unit(.5, "cm")
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::ggtitle(label)
}

grouped_sig_barplot <- function(hist_groups, sig_num_df, output_dir = getwd(),
                                label) {
  # Given the data.frame output from `calc_mut_per_sig` where the number of
  # mutations per Mb that belong to each sample x signature combo, plot the
  # number of mutations per Mb for each sample as a grouped bar plot for all
  # signatures.
  #
  # Args:
  #   hist_groups: a vector of `short_histology` groups to each be plotted
  #                individually.
  #   sig_num_df: a data.frame with number of mutations per Mb that belong to
  #               each sample x signature combo
  #   output_dir: where the plots should be saved. If this directory doesn' exist,
  #               will create it. Default is current directory.
  #   label: a character string for the title of the plot to be passed to it's
  #          png name and ggtitle.
  #
  # Returns:
  #  A grouped barplot saved as png with themutations per Mb for each sample
  #  from each signature. 

  # Make the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Do this for each histology group provided in. the vector
  for (hist_group in hist_groups) {
    # Narrow the df down to just this histology's group
    histology_df <- sig_num_df %>%
      dplyr::filter(short_histology == hist_group, mut_per_mb > 0) %>%
      dplyr::mutate("Signatures" = droplevels(variable))

    # Make the grouped bar plot
    histology_df %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(Tumor_Sample_Barcode, -mut_per_mb), y = mut_per_mb)) +
      ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = Signatures)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1)) +
      ggplot2::ylab("Mutations per Mb") +
      ggplot2::xlab("") +
      ggplot2::theme_classic() +
      ggplot2::ggtitle(paste(hist_group, label, "signatures"))

    # Save the plot
    ggplot2::ggsave(file.path(
      output_dir,
      paste0(
        "barplot_",
        hist_group,
        "_", label,
        "_mutation_sig.png"
      )
    ))
  }
}
