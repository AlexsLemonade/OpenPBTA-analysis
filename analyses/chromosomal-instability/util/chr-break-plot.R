# Functions for chromosomal instability calculations
#
# C. Savonen for ALSF - CCDL
#
# 2020

map_density_plot <- function(granges,
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
  #
  # For setting the scale later, need to get y's max
  max_y <- max(
    dplyr::pull(
      data.frame(granges@elementMetadata@listData),
      !!rlang::sym(y_val)
    )
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
    ggplot2::scale_y_continuous(breaks = seq(0, max_y, by = 2))

  # Print out plot
  density_plot@ggplot
}

chr_break_plot <- function(granges_list,
                           plot_name,
                           y_val) {
  # A wrapper function to make a 3 row chromosomal map plot for a set of GRanges
  # objects that contain common_density, cnv_density, and sv_density.
  #
  # Args:
  #   granges_list: A list of Granges object to plot as a combination plot
  #   plot_name: a character string specifying the plot
  #   y_val: to be passed to map_density plot for mapping. 
  #
  # Returns:
  #  ggplot of chromosomal mapping of the y value given.
  #
  # Make combined SV and CNV plot
  common_plot <- map_density_plot(granges_list$common_density,
    y_val = y_val,
    y_lab = "Breaks per Mb",
    color = "blue",
    main_title = "Common Breaks"
  )

  # Make CNV plot
  cnv_plot <- map_density_plot(granges_list$cnv_density,
    y_val = y_val,
    y_lab = "Breaks per Mb",
    color = "darkgreen",
    main_title = "CNV Breaks"
  )

  # Make SV plot
  sv_plot <- map_density_plot(granges_list$sv_density,
    y_val = y_val,
    y_lab = "Breaks per Mb",
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
    common_plot,
    cnv_plot,
    sv_plot,
    nrow = 4,
    axis = "b",
    rel_heights = c(2, 5, 5, 5)
  )

  # Save plot to PNG
  cowplot::save_plot(
    plot = full_plot,
    filename = file.path(plots_dir, paste0(plot_name, "_breaks.png")),
    base_height = 7,
    base_width = 20
  )

  # Print out the plot while we are here
  full_plot
}
