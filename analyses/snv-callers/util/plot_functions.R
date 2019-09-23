# Plotting functions for MAF data
#
# C. Savonen for ALSF - CCDL
# 2019
#
################################################################################
############################# Plotting Functions ###############################
################################################################################

base_change_plot <- function(vaf_df, exp_strategy = "BOTH", filter_cutoff = 0) {
  # Plot the number of base changes as a barplot. Need a MAF data.frame that
  # that has `base_change` column set up with `calculate_vaf.`
  #
  # Args:
  #   vaf_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `calculate_vaf`
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both` Case insensitive.
  #   filter_cutoff: A numeric number to only keep groups larger than it.
  #
  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy ifspecified
  if (exp_strategy != "BOTH") {
    vaf_df <- vaf_df %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Count the number of each type of base change
  base_count <- vaf_df %>%
    dplyr::group_by(change, experimental_strategy) %>%
    dplyr::summarise(base_count = dplyr::n()) %>%
    dplyr::filter(base_count > filter_cutoff)

  # Plot this as a barplot
  barplot <- ggplot2::ggplot(
    base_count,
    ggplot2::aes(x = reorder(change, -base_count), y = base_count)
  )

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    barplot <- barplot +
      ggplot2::geom_bar(
        position = "dodge", stat = "identity",
        fill = "navyblue"
      ) +
      ggplot2::theme(legend.position = "none")
  } else {
    barplot <- barplot +
      ggplot2::geom_bar(
        position = "dodge", stat = "identity",
        ggplot2::aes(fill = experimental_strategy)
      )
  }
  # Add some other aesthetics
  barplot <- barplot +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::xlab("") +
    ggplot2::ylab("Count") +
    colorblindr::scale_fill_OkabeIto()

  return(barplot)
}

depth_vs_vaf_plot <- function(vaf_df, exp_strategy = "BOTH") {
  # Plot read depth versus VAF
  #
  # Args:
  #   vaf_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `calculate_vaf`
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both` Case insensitive.
  #
  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy if specified
  if (exp_strategy != "BOTH") {
    vaf_df <- vaf_df %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Plot this as a scatterplot
  barplot <- ggplot2::ggplot(
    vaf_df,
    ggplot2::aes(x = log2(n_depth), y = vaf)
  ) +
    ggplot2::theme_classic()

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    barplot <- barplot +
      ggplot2::geom_point(color = "navyblue", alpha = 0.15) +
      ggplot2::theme(legend.position = "none")
  } else {
    barplot <- barplot +
      ggplot2::geom_point(ggplot2::aes(color = experimental_strategy), alpha = 0.15) +
      colorblindr::scale_color_OkabeIto()
  }
  return(barplot)
}

snv_region_plot <- function(maf_annot, exp_strategy = "BOTH", filter_cutoff = 0) {
  # Plot the genomic region analysis as a summarized barplot
  #
  # Args:
  #   maf_annot: MAF formatted data that has been turned into a data.frame and has
  #           been run through `annotr_maf`
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both`.
  #   filter_cutoff: A numeric number to only keep groups larger than it.
  #
  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy if specified
  if (exp_strategy != "BOTH") {
    maf_annot <- maf_annot %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Count the number of each type of base change
  type_count_df <- maf_annot %>%
    dplyr::group_by(type, experimental_strategy) %>%
    dplyr::summarise(type_count = dplyr::n()) %>%
    dplyr::filter(type_count > filter_cutoff)

  # Plot this as a barplot
  barplot <- ggplot2::ggplot(
    type_count_df,
    ggplot2::aes(x = reorder(type, -type_count), y = type_count)
  )

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    barplot <- barplot +
      ggplot2::geom_bar(
        position = "dodge", stat = "identity",
        fill = "navyblue"
      ) +
      ggplot2::theme(legend.position = "none")
  } else {
    barplot <- barplot +
      ggplot2::geom_bar(
        position = "dodge", stat = "identity",
        ggplot2::aes(fill = experimental_strategy)
      )
  }

  # Add some other aesthetics
  barplot <- barplot +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::xlab("") +
    ggplot2::ylab("Count") +
    colorblindr::scale_fill_OkabeIto()

  return(barplot)
}

cosmic_plot <- function(vaf_df, exp_strategy = "BOTH") {
  # Plot the VAF for COSMIC vs non-COSMIC mutations
  #
  # Args:
  #   vaf_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `calculate_vaf`
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both`.
  #
  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy if specified
  if (exp_strategy != "BOTH") {
    vaf_df <- vaf_df %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Get the overlap with COSMIC mutations using special function
  maf_cosmic <- suppressWarnings(find_cosmic_overlap(vaf_df))

  # Plot this as a violin plot
  barplot <- maf_cosmic %>%
    ggplot2::ggplot(ggplot2::aes(x = overlap_cosmic, y = vaf))

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    barplot <- barplot +
      ggplot2::geom_violin(fill = "navyblue") +
      ggplot2::theme(legend.position = "none")
  } else {
    barplot <- barplot +
      ggplot2::geom_violin(ggplot2::aes(fill = experimental_strategy))
  }
  # Add some aesthetics
  barplot <- barplot +
    ggplot2::theme_classic() +
    ggplot2::xlab("Same mutation and base change found in COSMIC") + 
    colorblindr::scale_fill_OkabeIto()

  return(barplot)
}

tmb_plot <- function(tmb_df, exp_strategy = "BOTH", x_axis = "short_histology") {
  # Plot Tumor Mutational Burden as a jitter plot by a variable of choice.
  #
  # Args:
  #   tmb_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `vaf`
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both`.
  #   x_axis: what variable you would like to use to plot on the x-axis
  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy ifspecified
  if (exp_strategy != "BOTH") {
    tmb_df <- tmb_df %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Plot this as a jitter plot
  barplot <- tmb_df %>%
    ggplot2::ggplot(ggplot2::aes(x = eval(parse(text = x_axis)), y = tmb))

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    barplot <- barplot +
      ggplot2::geom_jitter(color = "navyblue", alpha = 0.2) +
      ggplot2::theme(legend.position = "none")
  } else {
    barplot <- barplot +
      ggplot2::geom_violin(ggplot2::aes(fill = experimental_strategy))
  }

  # Add some aesthetics
  barplot <- barplot +
    ggplot2::theme_classic() +
    ggplot2::ylab("Tumor Mutational Burden") +
    ggplot2::xlab(x_axis) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + 
    colorblindr::scale_fill_OkabeIto()

  return(barplot)
}
