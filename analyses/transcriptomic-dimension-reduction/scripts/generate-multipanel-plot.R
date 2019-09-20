# Bethell and Taroni for CCDL 2019
# This script creates multipanel plots from plot lists saved as RDS files, the
# output of get-plot-list.R
#
# Command line usage:
# 
# 

library(ggplot2)

# read in plot
# plot_list <- readr::read_rds(plot_file)

plot_legend <- cowplot::get_legend(
  plot_list[[1]] +  theme(legend.box.margin = margin(0, 0, 0, 12))
) 

plot_list <- lapply(plot_list, 
                    function(p) p + theme(legend.position = "none"))

multipanel_plot <- cowplot::plot_grid(
  plotlist = plot_list,
  nrow = 1,
  align = "h",
  labels = "AUTO"
)