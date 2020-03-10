# Function for CDF plot
#
# C. Savonen for ALSF - CCDL
# 2020 

cdf_plot <- function(df,
                     plot_title,
                     num_col,
                     group_col,
                     color = "black",
                     n_group = 5,
                     x_lim = c(-1.2, 1.2),
                     y_lim,
                     x_lab,
                     y_lab, 
                     breaks) {
  # For a data.frame, create a Cummulative distribution function plot of the
  # `num_col` and `group_col` column data provided.
  # Groups in group_col will be plotted separately and reordered based on 
  # medians. Data is transformed to log10 for y-axis. 
  #
  # Args:
  #   df: a data.frame that includes the columns specified in `group_col` and
  #       `num_col`.
  #   plot_title: A string to be passed to ggplot2::ggtitle()
  #   group_col: a character string that indicates what column should be used to
  #              group and facet_wrap the CDF plot.
  #   num_col: a character string that indicates what column of numeric data
  #            should be plotted.
  #   color: What color should the points be. Will be passed to `geom_point`.
  #   n_group: What is the minimum number of samples needed when grouping by
  #            group_col to include in the plot? Default is 5.
  #   x_lim: a numeric vector of length two to be used as minimum and maximum x
  #          axis limits
  #   y_lim: a numeric vector of length two to be used as minimum and maximum y
  #          axis limits
  #   x_lab: a string to be passed to ggplot2::xlab() for an x-axis label
  #   y_lab: a string to be passed to ggplot2::ylab() for an y-axis label
  #   breaks: A numeric vector to be passed to ggplot2::scale_y_continuous for 
  #           marking ticks on the y axis. 
  #
  #
  # Returns:
  # A CDF plot of the num_col and group_col column data provided.

  # List the columns we need
  needed_cols <- c(num_col, group_col)

  # Get logical vector indicating which are in metadata
  found_cols <- (needed_cols %in% colnames(df))

  # If not all the columns are found, stop
  if (!all(found_cols)) {
    stop(cat(
      "The following column names specified for the CDF plot: \n",
      paste(needed_cols[which(found_cols)], collapse = "\n"),
      "\n ...were not found in the specified data.frame.",
      "Check your `num_col` and `group_col` arguments."
    ))
  }
  # Make sure this is a data.frame
  df <- df %>%
    as.data.frame(stringsAsFactors = FALSE)

  # Pull out this data
  group_var <- df %>%
    dplyr::pull(!!group_col)

  # Pull out this data
  num_var <- df %>%
    dplyr::pull(!!num_col)

  df <- df %>%
    # We only really need these two variables from data.frame
    dplyr::transmute(
      group = tools::toTitleCase(group_var),
      number = as.numeric(num_var)
    ) %>%
    # Group by specified column
    dplyr::group_by(group) %>%
    # Only keep groups with the specified minimum number of samples
    dplyr::filter(dplyr::n() > n_group) %>%
    # Calculate group median
    dplyr::mutate(
      group_median = median(number),
      group_rank = rank(number, ties.method = "first") / dplyr::n(),
      sample_size = paste0("n = ", dplyr::n())
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = reorder(group, group_median)) %>%
    # Now we will plot these as cumulative distribution plots
    ggplot2::ggplot(ggplot2::aes(
      x = group_rank,
      y = number
    )) +
    ggplot2::geom_point(color = color) +
    # Add summary line for median
    ggplot2::geom_segment(
      x = 0, xend = 1, color = "grey",
      ggplot2::aes(y = group_median, yend = group_median)
    ) +
    # Separate by histology
    ggplot2::facet_wrap(~ group + sample_size, nrow = 1, strip.position = "bottom") +
    ggplot2::theme_classic() +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    # Transform to log10 make non-log y-axis labels
    ggplot2::scale_y_continuous(
      trans = "log1p",
      limits = y_lim, 
      breaks = breaks
    ) +
    ggplot2::scale_x_continuous(limits = x_lim) +
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
