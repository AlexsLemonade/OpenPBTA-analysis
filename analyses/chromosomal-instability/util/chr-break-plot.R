# Functions for chromosomal instability calculations
#
# C. Savonen for ALSF - CCDL
#
# 2020

map_breaks_plot <- function(granges,
                            y_val,
                            y_lab,
                            color,
                            main_title) {
  # Given a GRanges object, plot it y value along its chromosomal mappings using
  # ggbio.
  #
  # Args:
  #   granges: A Granges object to plot
  #   y_val: a character string of the columnname in listData spot of the
  #          GenomicRanges to plot on the y axis
  #   color: a color parameter
  #   y_lab: a character string to use for the ylabel. Will be passed to
  #          ggplot2::ylab argument.
  #   main_title: a character string to use for the main title.
  #
  # Returns:
  #  ggplot of chromosomal mapping of the y value given.

  # For the y-axis ticks, default is to print every two
  by_interval <- 2

  # For setting the scale later, need to get y's max
  max_y <- max(
    data.frame(granges@elementMetadata@listData) %>%
      dplyr::pull(
        !!rlang::sym(y_val)
      ),
    na.rm = TRUE
  )
  # Make the density plot
  density_plot <- ggbio::autoplot(granges, ggplot2::aes(y = !!rlang::sym(y_val)),
    geom = "line", scales = "free_x", space = "free_x",
    colour = color
  ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 3, angle = 45, hjust = 1)) +
    ggplot2::ylab(y_lab) +
    ggplot2::ggtitle(main_title) +
    ggplot2::scale_y_continuous(breaks = seq(0, max_y, by = by_interval))

  # Print out plot
  density_plot@ggplot
}

multipanel_break_plot <- function(granges_list,
                                  plot_name,
                                  y_val,
                                  y_lab,
                                  plot_dir) {
  # A wrapper function to make a 3 row chromosomal map plot for a set of GRanges
  # objects that contain union_of_breaks, cnv_breaks, and sv_breaks.
  #
  # Args:
  #   granges_list: A list of Granges object to plot as a combination plot
  #   plot_name: a character string specifying the plot
  #   y_val: to be passed to map_breaks plot for mapping.
  #   y_lab: to be passed to map_breaks plot for y axis label
  #   plot_dir: a file path where you would like the plot PNG to be saved.
  #
  # Returns:
  #  ggplot of chromosomal mapping of the y value given.
  #
  # Make combined SV and CNV plot
  intersection_of_plot <- map_breaks_plot(granges_list$intersection_of_breaks,
    y_val = y_val,
    y_lab = y_lab,
    color = "blue",
    main_title = "Intersection of Breaks"
  )

  # Make CNV plot
  cnv_plot <- map_breaks_plot(granges_list$cnv_breaks,
    y_val = y_val,
    y_lab = y_lab,
    color = "darkgreen",
    main_title = "CNV Breaks"
  )

  # Make SV plot
  sv_plot <- map_breaks_plot(granges_list$sv_breaks,
    y_val = y_val,
    y_lab = y_lab,
    color = "orange",
    main_title = "SV Breaks"
  )

  # Make a title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(paste(plot_name, " - Chromosomal Break Density"),
      fontface = "bold", x = .4, hjust = 0, size = 12
    )

  # Put all plots and title together
  full_plot <- cowplot::plot_grid(title,
    intersection_of_plot,
    cnv_plot,
    sv_plot,
    nrow = 4,
    axis = "b",
    rel_heights = c(2, 5, 5, 5)
  )

  # Save plot to PNG
  cowplot::save_plot(
    plot = full_plot,
    filename = file.path(plot_dir, paste0(plot_name, "_breaks.png")),
    base_height = 7,
    base_width = 20
  )
}

breaks_cdf_plot <- function(density_file, metadata_df, cdf_plot_file) {
  # Given a genome wide breaks density file path, plot the CDF distribution 
  # for it by histology.
  #
  # Args:
  #   density_file: a chromosomal breaks density file path where each sample is
  #                 a row with these columns: samples, experimental_strategy,
  #                 breaks_count. 
  #   metadata_df: a data.frame with Kids_First_Biospecimen_ID, short_histology
  #             columns
  #   cdf_plot_file: A file.path to where the plot should be saved
  #   
  # Returns:
  #  Saves the CDF breaks plot to a file and also prints it out. 
  
  # Read in the breaks file
  breaks_cdf <- readr::read_tsv(density_file) %>%
    # Add on the short histology column from the metadata for plotting purposes
    dplyr::inner_join(metadata_df %>%
                        dplyr::select(Kids_First_Biospecimen_ID, short_histology),
                      by = c("samples" = "Kids_First_Biospecimen_ID")
    ) %>%
    # Only plot histologies groups with more than `min_samples` number of samples
    dplyr::group_by(tools::toTitleCase(short_histology), add = TRUE) %>%
    # Only keep groups with this amount of samples
    dplyr::filter(dplyr::n() > params$min_samples) %>%
    # Calculate histology group mean
    dplyr::mutate(
      hist_mean = mean(breaks_count),
      hist_rank = rank(breaks_count, ties.method = "first") / dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    # Now we will plot these as cummulative distribution plots
    ggplot2::ggplot(ggplot2::aes(
      x = hist_rank, 
      y = breaks_count,
      color = short_histology
    )) +
    ggplot2::geom_point() + 
    # Add summary line for mean
    ggplot2::geom_segment(x = 0, xend = 1, color = "grey", 
                          ggplot2::aes(y = hist_mean, yend = hist_mean)) +
    # Separate by histology
    ggplot2::facet_wrap(~ short_histology, nrow = 1, strip.position = "bottom") +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("Chromosomal Breaks per Genome") +
    # Transform to log10 make non-log y-axis labels
    ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30)) +
    ggplot2::scale_x_continuous(limits = c(-0.2, 1.2), breaks = c()) +
    # Making it pretty
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 10, angle = 90, hjust = 1), 
      strip.background = ggplot2::element_rect(fill = NA, color = NA)
    )
  
  # Save the plot to a png
  ggplot2::ggsave(cdf_plot_file,
                  plot = breaks_cdf, width = 40, height = 20, unit = "cm")
  
  return(breaks_cdf)
}
