# Functions for chromosomal instability calculations
#
# C. Savonen for ALSF - CCDL
#
# 2020

sina_plot <- function(data, 
                      colour = "#630882", 
                      x_axis,  
                      y_axis, 
                      title, 
                      y_label) {
  # Given a data.frame with metadata and `tmb` information, make a sina plot from it
  # plots log10 of the data. 
  #
  # Args:
  #
  # data : a data.frame with `tmb` as columns
  # x_axis : A name of a column in data that you would like to make the x_axis from
  #           Given as a character string. 
  # y_axis : A name of a column in data that you would like to make the x_axis from
  #           Given as a character string. 
  # colour: the color you would like the points to be
  # title: The title for the plot. To be passed to `ggtitle`
  # y_label: Label for y-axis. To be passed to `ggplot2::ylab`
  #
  # Returns:
  #   A sina plot by histology group
  data %>%
    # eval parse bit will treat the variable as the text it has stored
    ggplot2::ggplot(ggplot2::aes(x = reorder(!!rlang::sym(x_axis), !!rlang::sym(y_axis), mean), 
                                 y = log10(!!rlang::sym(y_axis)))) +
    ggforce::geom_sina(color = colour, alpha = 0.3, maxwidth = 0.9) +
    ggplot2::stat_summary(
      fun.y = mean, fun.ymin = mean, fun.ymax = mean,
      geom = "crossbar",
      width = 0.95, size = 0.15
    ) +
    ggplot2::theme_classic() +
    ggplot2::ylab(y_label) +
    ggplot2::xlab("") +
    ggplot2::ylim(c(-3, 3)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::ggtitle(title) + 
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(colour = c("gray93", "white"), size = 9)
    )
}
