# Bethell and Taroni for CCDL 2019
# This script creates multipanel plots from plot lists saved as RDS files, the
# output of get-plot-list.R.
#
# Command line usage:
<<<<<<< HEAD
# 
# Rscript scripts/generate-multipanel-plot.R \
#   --plot_rds plots/plot_data/rsem_all_broad_histology_multiplot_list.RDS \
#   --plot_directory plots 
# 
=======
#
# Rscript scripts/generate-multipanel-plot.R \
#   --plot_rds plots/plot_data/rsem_all_broad_histology_multiplot_list.RDS \
#   --plot_directory plots
#
>>>>>>> upstream/master

# we're going to use this quite a bit
library(ggplot2)

#### command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--plot_rds"),
    type = "character",
    default = NULL,
    help = "RDS file that contains the plot list",
  ),
  optparse::make_option(
    c("-d", "--plot_directory"),
    type = "character",
    default = NULL,
    help = "output directory for plots"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

plot_rds <- opt$plot_rds
plot_directory <- opt$plot_directory
if (!dir.exists(plot_directory)) {
  dir.create(plot_directory, recursive = TRUE)
}

#### Make the multipanel plot --------------------------------------------------

<<<<<<< HEAD
# First, we're going to infer some things about the plot from the name of 
=======
# First, we're going to infer some things about the plot from the name of
>>>>>>> upstream/master
# the input file -- specifically, the plot title and the output file name

# Remove the path information + change file extension
rds_filename <- sub(".*/", "", plot_rds)
output_file <- sub("_multiplot_list.RDS", ".pdf", rds_filename)

# extract the method and RNA_library information -- we'll use this as a title
quant_method <- stringr::word(output_file, 1, sep = "_")
RNA_library <- stringr::word(output_file, 2, sep = "_")
<<<<<<< HEAD

# putting the title together
plot_title <- toupper(paste(quant_method, RNA_library))
=======
transform_method <- stringr::word(output_file, 3, sep = "_")

# putting the title together
plot_title <- toupper(paste(quant_method, RNA_library, transform_method))
>>>>>>> upstream/master

# We need to read in the list generated via get-plot-list.R
plot_list <- readr::read_rds(plot_rds)

# Extract legend from the plot list + adjust text size
plot_legend <- cowplot::get_legend(
  plot_list[[1]] +
    guides(color = guide_legend(ncol = 1)) +
    theme(text = element_text(size = 10),
          legend.box.margin = margin(15, 15, 15, 15))
<<<<<<< HEAD
) 

# We'll remove the legends from all the plots for arranging in a grid
plot_list <- lapply(plot_list, 
                    function(p) p + theme(legend.position = "none"))

# Multipanel plot of the dimension reduction scatter plots with panel labels
multipanel_plot <- cowplot::plot_grid(plotlist = plot_list, nrow = 1, 
                                       align = "h", labels = "AUTO")

# This is the plot title
plot_title_panel <- cowplot::ggdraw() + 
  cowplot::draw_label(plot_title, fontface = "bold", size = 15)

# Add the title such that it is centered over the scatter plots
multipanel_with_title <- cowplot::plot_grid(plot_title_panel, multipanel_plot, 
                                            ncol = 1, rel_heights = c(0.1, 1))

# Pad the right side of the multipanel plot that now has a title
multipanel_spacer <- cowplot::plot_grid(multipanel_with_title, NULL, 
                                        nrow = 1, rel_widths = c(1, 0.05))

# Add in the legend to the right of the multipanel plot
final_plot <- cowplot::plot_grid(multipanel_with_title, plot_legend, 
=======
)

# We'll remove the legends from all the plots for arranging in a grid
plot_list <- lapply(plot_list,
                    function(p) p + theme(legend.position = "none"))

# Multipanel plot of the dimension reduction scatter plots with panel labels
multipanel_plot <- cowplot::plot_grid(plotlist = plot_list, nrow = 1,
                                       align = "h", labels = "AUTO")

# This is the plot title
plot_title_panel <- cowplot::ggdraw() +
  cowplot::draw_label(plot_title, fontface = "bold", size = 15)

# Add the title such that it is centered over the scatter plots
multipanel_with_title <- cowplot::plot_grid(plot_title_panel, multipanel_plot,
                                            ncol = 1, rel_heights = c(0.1, 1))

# Pad the right side of the multipanel plot that now has a title
multipanel_spacer <- cowplot::plot_grid(multipanel_with_title, NULL,
                                        nrow = 1, rel_widths = c(1, 0.05))

# Add in the legend to the right of the multipanel plot
final_plot <- cowplot::plot_grid(multipanel_with_title, plot_legend,
>>>>>>> upstream/master
                                 nrow = 1, rel_widths = c(1, 0.25)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# save to file
<<<<<<< HEAD
cowplot::save_plot(file.path(plot_directory, output_file), final_plot, 
=======
cowplot::save_plot(file.path(plot_directory, output_file), final_plot,
>>>>>>> upstream/master
                   base_width = 21, base_height = 7)
