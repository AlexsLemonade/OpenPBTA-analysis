# Bethell and Taroni for CCDL 2019
# This script generates a list of scatter plots, saved as an RDS. This plot
# list can be used downstream for multigrid plots.
# 
# Command line usage:
#
# Rscript --vanilla scripts/get-plot-list.R \
#   --input_directory results \
#   --filename_lead rsem_all \
#   --output_directory plots \
#   --color_variable RNA_library \

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-i", "--input_directory"),
    type = "character",
    default = NULL,
    help = "directory that contains input fiiles",
  ),
  optparse::make_option(
    c("-f", "--filename_lead"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to identify input files"
  ),
  optparse::make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "output directory"
  ),
  optparse::make_option(
    c("-c", "--color_variable"),
    type = "character",
    default = NULL,
    help = "Name of column that will be used to color points in plots"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

input_directory <- opt$input_directory
filename_lead <- opt$filename_lead
color_variable <- opt$color_variable

# Create the output directory if it does not exist
output_directory <- opt$output_directory 
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

#### Function ------------------------------------------------------------------

# Detect the ".git" folder -- this will be in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# source the functions for dimension reduction
source(file.path(root_dir, "analyses", "transcriptomic-dimension-reduction", 
                 "util", "dimension-reduction-functions.R"))

#### Files ---------------------------------------------------------------------

pca_file <- file.path(input_directory,
                      paste0(filename_lead, "_pca_scores_aligned.tsv"))
umap_file <- file.path(input_directory,
                       paste0(filename_lead, "_umap_scores_aligned.tsv"))
tsne_file <- file.path(input_directory,
                       paste0(filename_lead, "_tsne_scores_aligned.tsv"))

# this will be a logical that tells us whether or not t-SNE was run via
# run-dimension-reduction.R script
tsne_run <- file.exists(tsne_file)

# the output of this will be an RDS that contains a list of plots
output_file <- file.path(output_directory,
                         paste(filename_lead, color_variable,
                               "multiplot_list.RDS", sep = "_"))

#### Generate plot lists -------------------------------------------------------

# empty list that will hold plots
plot_list <- list()

pca_df <- readr::read_tsv(pca_file)
plot_list[["PCA"]] <- plot_dimension_reduction(aligned_scores_df = pca_df,
                                               point_color = color_variable,
                                               x_label = "PC1",
                                               y_label = "PC2")

umap_df <- readr::read_tsv(umap_file)
plot_list[["UMAP"]] <- plot_dimension_reduction(aligned_scores_df = umap_df,
                                                point_color = color_variable,
                                                x_label = "UMAP1",
                                                y_label = "UMAP2")

# if t-SNE was run, add that plot to the list
if (tsne_run) {
  tsne_df <- readr::read_tsv(tsne_file)
  plot_list[["t-SNE"]] <- plot_dimension_reduction(aligned_scores_df = tsne_df,
                                                   point_color = color_variable,
                                                   x_label = "t-SNE1",
                                                   y_label = "t-SNE2")
}

# save plot list as an RDS
readr::write_rds(plot_list, path = output_file)
