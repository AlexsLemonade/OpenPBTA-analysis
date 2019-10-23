# 01-generate-independent-specimens.R


# Base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "independent-samples")

# Load the optparse library
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# source sample selection function
source(file.path(analysis_dir, "independent-samples.R"))

set.seed(201910)

# Parse options

option_list <- list(
  make_option(
    c("-f", "--histology_file"),
    type = "character",
    default = NULL,
    help = "path to the histology tsv file, relative to project root",
  ),
  make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "path to output directory, relative to project root"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# set output files


# Read histology file
sample_df <- readr::read_tsv(file.path(root_dir, opts$histology_file))


# Filter to only WGS samples from tumors, where composition is known to be Solid Tissue
# Note that there are some samples with unknown composition, but these will be ignored for now.
wgs_samples <- sample_df %>%
  dplyr::filter(experimental_strategy == "WGS",
                sample_type == "Tumor", 
                composition == "Solid Tissue")

primary_only <- independent_samples(wgs_samples, tumor_types = "primary")

primary_plus <- independent_samples(wgs_samples, tumor_types = "prefer_primary")


         
  
