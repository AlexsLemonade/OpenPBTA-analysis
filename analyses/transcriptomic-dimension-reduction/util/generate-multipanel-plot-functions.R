# Intended for import only
# Chante Bethell and Jaclyn Taroni for CCDL 2019
#
# Function to generate a multipanel plot
generate_multipanel_plot <- function(plot_list,
                                     plot_title,
                                     output_directory,
                                     output_filename) {
  # Given two plots generated in a previous function, make a multipanel plot,
  # save, and display.
  #
  # Args:
  #   plot_list: list of individual plots to be put together
  #   plot_title: string to add as title of output plot
  #   output_directory: the path to the directory where the multipanel plot
  #                     should be saved
  #   output_filename: the filename for the output multipanel plot
  #
  # Returns:
  #   final_plot: multipanel plot with visualizations of the individual plots
  
  # Extract legend from the plot list + adjust text size
  plot_legend <- cowplot::get_legend(
    plot_list[[1]] +
      guides(color = guide_legend(ncol = 1)) +
      theme(
        text = element_text(size = 10),
        legend.box.margin = margin(15, 15, 15, 15)
      )
  )
  
  # We'll remove the legends from all the plots for arranging in a grid
  plot_list <- lapply(plot_list,
                      function(p)
                        p + theme(legend.position = "none"))
  
  # Multipanel plot of the dimension reduction scatter plots with panel labels
  multipanel_plot <-
    cowplot::plot_grid(
      plotlist = plot_list,
      nrow = 1,
      align = "h",
      labels = "AUTO"
    )
  
  # This is the plot title
  plot_title_panel <- cowplot::ggdraw() +
    cowplot::draw_label(plot_title, fontface = "bold", size = 15)
  
  # Add the title such that it is centered over the scatter plots
  multipanel_with_title <-
    cowplot::plot_grid(
      plot_title_panel,
      multipanel_plot,
      ncol = 1,
      rel_heights = c(0.1, 1)
    )
  
  # Pad the right side of the multipanel plot that now has a title
  multipanel_spacer <- cowplot::plot_grid(multipanel_with_title,
                                          NULL,
                                          nrow = 1,
                                          rel_widths = c(1, 0.05))
  
  # Add in the legend to the right of the multipanel plot
  final_plot <- cowplot::plot_grid(
    multipanel_with_title,
    plot_legend,
    nrow = 1,
    rel_widths = c(1, 0.25)
  ) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  # Save to file
  cowplot::save_plot(
    file.path(output_directory, output_filename),
    final_plot,
    base_width = 24,
    base_height = 8
  )
  
  # Display final plot
  final_plot
}
