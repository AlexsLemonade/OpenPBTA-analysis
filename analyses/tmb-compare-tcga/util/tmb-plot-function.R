# Function for TMB CDF plot
#
# C. Savonen for ALSF - CCDL

tmb_cdf_plot <- function(tmb_df,
                         plot_title,
                         color, 
                         n_group) {
  # Print CDF TMB plot
  #
  # Args:
  #   tmb_df: a data.frame with `tmb` and `short_histology` as column names.
  #   plot_title: A string to be passed to "ggtitle".
  #   color: What color should the points be. Will be passed to `geom_point`. 
  #   n_group: What is the minimum number of samples needed in a histology group 
  #            for it to be plotted. 
  #
  # Returns:
  # A CDF plot of TMB with histology
  tmb_df %>%
    as.data.frame() %>%
    dplyr::mutate(short_histology = tools::toTitleCase(short_histology)) %>%
    # Only plot histologies groups with more than `min_samples` number of samples
    dplyr::group_by(short_histology, add = TRUE) %>%
    # Only keep groups with this amount of samples
    dplyr::filter(dplyr::n() > n_group) %>%
    # Calculate histology group mean
    dplyr::mutate(
      hist_median = median(tmb),
      hist_rank = rank(tmb, ties.method = "first") / dplyr::n(),
      sample_size = paste0("n = ", dplyr::n())
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(short_histology = reorder(short_histology, hist_median)) %>%
    # Now we will plot these as cummulative distribution plots
    ggplot2::ggplot(ggplot2::aes(
      x = hist_rank,
      y = tmb
    )) +
    ggplot2::geom_point(color = color) +
    # Add summary line for mean
    ggplot2::geom_segment(
      x = 0, xend = 1, color = "grey",
      ggplot2::aes(y = hist_median, yend = hist_median)
    ) +
    # Separate by histology
    ggplot2::facet_wrap(~ short_histology + sample_size, nrow = 1, strip.position = "bottom") +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("Coding mutations per Mb") +
    # Transform to log10 make non-log y-axis labels
    ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 1000, 10000, 30000),
                                limits = c(0, 30000)) +
    ggplot2::scale_x_continuous(limits = c(-0.2, 1.2), breaks = c()) +
    # Making it pretty
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 10, angle = 90, hjust = 1),
      strip.background = ggplot2::element_rect(fill = NA, color = NA)
    ) +
    ggplot2::ggtitle(plot_title)
}
