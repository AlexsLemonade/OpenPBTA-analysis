# Take in oncoprint-goi-lists-OpenPBTA.csv and create a goi file for each
# column with associated genes of interest for each specified broad histology

# Use a wildcard to find the list for the appropriate histology in shell script or `01-plot.R`

# Chante Bethell for CCDL 2021
#
# # #### USAGE
# This script is intended to be sourced in the script as follows:
# 
# Rscript --vanilla util/prepare-goi-lists.R \
#   --goi_file data/oncoprint-goi-lists-OpenPBTA.csv \
#   --broad_histology "Low-grade astrocytic tumor"

#### Set Up --------------------------------------------------------------------

library(dplyr)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to output directory for plots produced
data_dir <-
  file.path(root_dir, "analyses", "oncoprint-landscape", "data")

if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-g", "--goi_file"),
    type = "character",
    default = NULL,
    help = "file path to comma-separatedGOI file that contains genes of interest per broad histology",
  ),
  optparse::make_option(
    c("-b", "--broad_histology"),
    type = "character",
    default = NULL,
    help = "name of `broad_histology` value to prepare associated goi list for"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

all_goi_df <- readr::read_csv(opt$goi_file)

prepare_histology_list <- function(goi_df,
                                   broad_histology) {
  # Given a data.frame with columns of histology specific genes of interest
  # and a specific broad histology, return a TSV file with the genes of
  # interest for the specified broad histology.
  #
  # Args:
  #   goi_df: data.frame with genes of interest separated by broad histology
  #           and supported by sources
  #   broad_histology: broad histology string whose list we want to prepare
  
  # Match the desired broad histology label with the short histology column label
  if (broad_histology == "Low-grade astrocytic tumor") {
    short_histology <- "LGAT"
  }
  else if (broad_histology == "Diffuse astrocytic and oligodendroglial tumor") {
    short_histology <- "HGAT"
  }
  else if (broad_histology == "Other CNS") {
    short_histology <- "Other"
  }
  else {
    short_histology = broad_histology
  }
  
  # Filter the data.frame and write a TSV file containing the specified broad
  # histology genes of interest
  histology_goi_list <- all_goi_df %>%
    select(short_histology) %>%
    tidyr::drop_na() %>%
    unique() %>%
    readr::write_tsv(file.path(data_dir, paste0(tolower(
      gsub(" ", "-", broad_histology)
    ), "_goi_list.tsv")))
}

# Run the custom `prepare_histology_list` function
histology_list <- prepare_histology_list(goi_df = all_goi_df,
                                         broad_histology = opt$broad_histology)
