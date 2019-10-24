# 01-generate-independent-specimens.R
#
# Josh Shapiro for CCDL 2019
# 
# Purpose: Generate tables of independent specimens where no two specimens are 
#   chosen from the same individual.
#
# Option descriptions
# -f, --histology_file : File path to where you would like the annotation_rds file to be
#               stored
# -o,--output_directory : Output directory
# 
# example invocation:
# Rscript analyses/independent-samples/01-generate-independent-specimens.R \
#   -f data/pbta-histologies.tsv \
#   -o analyses/independent-samples/results


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
out_dir <- file.path(root_dir, opts$output_directory)
if (!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE)
}

wgs_primary_file <- file.path(out_dir, 
                              "independent-specimens.wgs.primary.tsv")
wgs_primplus_file <- file.path(out_dir, 
                               "independent-specimens.wgs.primary-plus.tsv")
wgswxs_primary_file <- file.path(out_dir, 
                                 "independent-specimens.wgswxs.primary.tsv")
wgswxs_primplus_file <- file.path(out_dir, 
                                  "independent-specimens.wgswxs.primary-plus.tsv")

# Read histology file
sample_df <- readr::read_tsv(file.path(root_dir, opts$histology_file), 
                             col_types = readr::cols()) # suppress parse message


# Filter to only samples from tumors, where composition is known to be Solid Tissue
# Note that there are some samples with unknown composition, but these will be ignored for now.
tumor_samples <- sample_df %>%
  dplyr::filter(sample_type == "Tumor", 
                composition == "Solid Tissue", 
                experimental_strategy %in% c("WGS", "WXS"))

# Generate WGS independent samples
wgs_samples <- tumor_samples %>%
  dplyr::filter(experimental_strategy == "WGS")

wgs_primary <- independent_samples(wgs_samples, tumor_types = "primary")
wgs_primary_plus <- independent_samples(wgs_samples, tumor_types = "prefer_primary")

# Generate lists for WXS only samples
wxs_only_samples <-  tumor_samples %>% 
  dplyr::filter(!(Kids_First_Participant_ID %in% 
                    wgs_samples$Kids_First_Participant_ID))

wxs_primary <- independent_samples(wxs_only_samples, tumor_types = "primary")
wxs_primary_plus <- independent_samples(wxs_only_samples, tumor_types = "prefer_primary")

# write files
message(paste(nrow(wgs_primary), "WGS primary specimens"))
readr::write_tsv(wgs_primary, wgs_primary_file)

message(paste(nrow(wgs_primary_plus), "WGS specimens (including non-primary)"))
readr::write_tsv(wgs_primary_plus, wgs_primplus_file)

message(paste(nrow(wgs_primary) + nrow(wxs_primary), "WGS+WXS primary specimens"))
readr::write_tsv(dplyr::bind_rows(wgs_primary, wxs_primary),
                 wgswxs_primary_file)

message(paste(nrow(wgs_primary_plus) + nrow(wxs_primary_plus), "WGS+WXS specimens (including non-primary)"))
readr::write_tsv(dplyr::bind_rows(wgs_primary_plus, wxs_primary_plus),
                 wgswxs_primplus_file)
