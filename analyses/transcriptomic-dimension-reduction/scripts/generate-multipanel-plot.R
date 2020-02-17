# Bethell and Taroni for CCDL 2019
# This script creates multipanel plots from plot lists saved as RDS files, the
# output of get-plot-list.R.
#
# Command line usage:
#
# Rscript scripts/generate-multipanel-plot.R \
#   --plot_rds plots/plot_data/rsem_all_broad_histology_multiplot_list.RDS \
#   --plot_directory plots
#

# we're going to use this quite a bit
library(ggplot2)

#### Command line options ------------------------------------------------------

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
  ),
  optparse::make_option(
    c("-f", "--function_script"),
    type = "character",
    default = NULL,
    help = "R script that contains custom function to generate multipanel plot"
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

#### Source custom function ---------------------------------------------------

# Source script that contains the custom `generate-multipanel-plot` function
source(opt$function_script)

#### Make the multipanel plot --------------------------------------------------

# First, we're going to infer some things about the plot from the name of
# the input file -- specifically, the plot title and the output file name

# Remove the path information + change file extension
rds_filename <- sub(".*/", "", plot_rds)
output_file <- sub("_multiplot_list.RDS", ".pdf", rds_filename)

# Extract the method and RNA_library information -- we'll use this as a title
quant_method <- stringr::word(output_file, 1, sep = "_")
RNA_library <- stringr::word(output_file, 2, sep = "_")
transform_method <- stringr::word(output_file, 3, sep = "_")

# putting the title together
plot_title <- toupper(paste(quant_method, RNA_library, transform_method))

# We need to read in the list of plots
plot_list <- readr::read_rds(plot_rds)

# Run the `generate-multipanel-plot` function
generate_multipanel_plot(plot_list = plot_list,
                         plot_title = plot_title,
                         output_directory = plot_directory,
                         output_filename = output_file)
