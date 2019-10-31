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
  # that has `vaf`, `base_change` and `experimental_strategy` columns that are
  # added with the `set_up_maf` function.
  #
  # Args:
  #   vaf_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `set_up_maf`.
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #                 `both` Case insensitive.
  #   filter_cutoff: A numeric number to only keep groups larger than it.
  #
  # Returns:
  # A barplot with the number of mutations with each type of base change noted
  # in the `base_change` column made in `set_up_maf`

  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy if specified
  if (exp_strategy != "BOTH") {
    vaf_df <- vaf_df %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Count the number of each type of base change
  base_count_df <- vaf_df %>%

    dplyr::mutate(
      change = as.factor(change),
      # Change factor level order so ins and del are at the end
      change = forcats::fct_relevel(change, "ins", "del", "long_change", after = Inf)
    ) %>%
    dplyr::group_by(change, experimental_strategy) %>%
    dplyr::summarise(base_count = dplyr::n()) %>%
    dplyr::filter(base_count > filter_cutoff)

  # Plot this as a barplot
  barplot <- ggplot2::ggplot(
    base_count_df,
    ggplot2::aes(
      x = change,
      y = base_count
    )
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
    ggplot2::xlab("Base Change") +
    ggplot2::ylab("Count") +
    colorblindr::scale_fill_OkabeIto()

  return(barplot)
}

depth_vs_vaf_plot <- function(vaf_df, exp_strategy = "BOTH") {
  # Plot read depth versus variant allele fraction. If both experimental strategies
  # are included, the plot is made color-coded by strategy.
  #
  # Args:
  #   vaf_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `set_up_maf` with metadata (specifically the
  #           `experimental_strategy` column added.
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
  scatterplot <- ggplot2::ggplot(
    vaf_df,
    ggplot2::aes(x = log2(n_depth), y = vaf)
  ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Read Depth")

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    scatterplot <- scatterplot +
      ggplot2::geom_point(color = "navyblue", alpha = 0.15) +
      ggplot2::theme(legend.position = "none")
  } else {
    scatterplot <- scatterplot +
      ggplot2::geom_point(ggplot2::aes(color = experimental_strategy), alpha = 0.15) +
      colorblindr::scale_color_OkabeIto()
  }
  return(scatterplot)
}

snv_region_plot <- function(maf_annot, exp_strategy = "BOTH", filter_cutoff = 0) {
  # Plot the genomic region analysis from the annotation added from `annotr_maf`
  # as a barplot with the number of mutations found in each type of genomic region.
  #
  # Args:
  #   maf_annot: MAF formatted data that has been turned into a data.frame and has
  #           been run through `annotr_maf`
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both`.
  #   filter_cutoff: A numeric number to only keep groups larger than it.
  #
  # Returns:
  # A barplot with the number of mutations that are found within each type of
  # genomic region in the `type` column
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
    ggplot2::aes(x = type, y = type_count)
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

cosmic_plot <- function(vaf_df, exp_strategy = "BOTH", cosmic_clean_file = NULL) {
  # Plot the VAF for COSMIC vs non-COSMIC mutations
  #
  # Args:
  #   vaf_df: Need a MAF data.frame that has `vaf`, `base_change` and
  #           `experimental_strategy` columns that are added with the
  #            `set_up_maf` function.
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both`.
  #   cosmic_clean_file: a file path to a TSV of COSMIC mutations that has been
  #                     been cleaned up to have the genomic coordinates separated
  #                     into Chr, Start, and End columns. This is passed to the
  #                     `find_cosmic_overlap` function.
  #
  # Returns:
  #  A violin plot with VAF plotted by whether or not it is overlapping with the
  #  COSMIC mutation set.
  #
  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy if specified
  if (exp_strategy != "BOTH") {
    vaf_df <- vaf_df %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Get the overlap with COSMIC mutations using special function
  maf_cosmic <- find_cosmic_overlap(vaf_df, cosmic_clean_file)

  # Plot this as a violin plot
  vioplot <- maf_cosmic %>%
    ggplot2::ggplot(ggplot2::aes(x = overlap_cosmic, y = vaf))

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    vioplot <- vioplot +
      ggplot2::geom_violin(fill = "navyblue") +
      ggplot2::theme(legend.position = "none")
  } else {
    vioplot <- vioplot +
      ggplot2::geom_violin(ggplot2::aes(fill = experimental_strategy))
  }
  # Add some aesthetics
  vioplot <- vioplot +
    ggplot2::theme_classic() +
    ggplot2::xlab("Overlaps with COSMIC mutation") +
    colorblindr::scale_fill_OkabeIto()

  return(vioplot)
}

tmb_plot <- function(tmb_df, exp_strategy = "BOTH", x_axis = "short_histology") {
  # Plot Tumor Mutational Burden as a jitter plot by a variable of choice.
  #
  # Args:
  #   tmb_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `vaf`
  #   exp_strategy: argument to specify whether to plot `wgs`, `wxs` samples or
  #   `both`.
  #   x_axis: what variable you would like to use to plot on the x-axis. Default
  #           is `short_histology` column.
  #
  # Returns:
  # A jitterplot that plots the TMB stats by the argument specified in x_axis
  # argument.
  #
  # Make this argument case insensitive
  exp_strategy <- toupper(exp_strategy)

  # Filter out based on experimental stategy ifspecified
  if (exp_strategy != "BOTH") {
    tmb_df <- tmb_df %>%
      dplyr::filter(experimental_strategy == exp_strategy)
  }

  # Plot this as a jitter plot
  jitterplot <- tmb_df %>%
    ggplot2::ggplot(ggplot2::aes(x = eval(parse(text = x_axis)), y = tmb))

  # Get rid of legend if both data aren't being plotted
  if (exp_strategy != "BOTH") {
    jitterplot <- jitterplot +
      ggplot2::geom_jitter(color = "navyblue", alpha = 0.2) +
      ggplot2::theme(legend.position = "none")
  } else {
    jitterplot <- jitterplot +
      ggplot2::geom_jitter(ggplot2::aes(color = experimental_strategy))
  }

  # Add some aesthetics
  jitterplot <- jitterplot +
    ggplot2::theme_classic() +
    ggplot2::ylab("Tumor Mutational Burden") +
    ggplot2::xlab(x_axis) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    colorblindr::scale_color_OkabeIto()

  return(jitterplot)
}
