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

# Load the libraries
library(magrittr)
library(optparse)


# source sample selection function
source(file.path(analysis_dir, "independent-samples.R"))

set.seed(201910)

# Parse options
option_list <- list(
  make_option(
    c("-f", "--histology_file"),
    type = "character",
    default = NULL,
    help = "path to the histology tsv file",
  ),
  make_option(
    c("-o", "--output_directory"),
    type = "character",
    default = NULL,
    help = "path to output directory"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# set output files
out_dir <- opts$output_directory
if (!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE)
}

wgs_primary_file <- file.path(out_dir, 
                              "independent-specimens.wgs.primary.tsv")
wgs_secondary_file <- file.path(out_dir, 
                              "independent-specimens.wgs.secondary.tsv")
wgs_primplus_file <- file.path(out_dir, 
                              "independent-specimens.wgs.primary-plus.tsv")
wgswxspanel_primary_file <- file.path(out_dir, 
                              "independent-specimens.wgswxspanel.primary.tsv")
wgswxspanel_secondary_file <- file.path(out_dir, 
                              "independent-specimens.wgswxspanel.secondary.tsv")
wgswxspanel_primplus_file <- file.path(out_dir, 
                              "independent-specimens.wgswxspanel.primary-plus.tsv")

# Read histology file
sample_df <- readr::read_tsv(opts$histology_file, 
                             guess_max = 100000,
                             col_types = readr::cols()) # suppress parse message


# Filter to only samples from tumors, where composition is known to be Solid Tissue
# Note that there are some samples with unknown composition, but these will be ignored for now.
tumor_samples <- sample_df %>%
  dplyr::filter(sample_type == "Tumor", 
                composition == "Solid Tissue", 
                experimental_strategy %in% c("WGS", "WXS", "Targeted Sequencing", "Targeted-Capture"))

# Generate WGS independent samples
wgs_samples <- tumor_samples %>%
  dplyr::filter(experimental_strategy == "WGS")

wgs_primary <- independent_samples(wgs_samples, tumor_types = "primary", seed = 2020)
wgs_secondary <- independent_samples(wgs_samples, tumor_types = "secondary", seed = 2020)
wgs_primary_plus <- independent_samples(wgs_samples, tumor_types = "prefer_primary", seed = 2020)

# Generate lists for WXS and Panel samples 
# WGS is generally preferred, so we will only include those where WGS is not available
wxs_panel_samples <-  tumor_samples %>% 
  dplyr::filter(!(Kids_First_Participant_ID %in% 
                  wgs_samples$Kids_First_Participant_ID))

wxs_panel_primary <- independent_samples(wxs_panel_samples, tumor_types = "primary", seed = 2020)
wxs_panel_secondary <- independent_samples(wxs_panel_samples, tumor_types = "secondary", seed = 2020)
wxs_panel_primary_plus <- independent_samples(wxs_panel_samples, tumor_types = "prefer_primary", seed = 2020)

# write files
message(paste(nrow(wgs_primary), "WGS primary specimens"))
readr::write_tsv(wgs_primary, wgs_primary_file)

message(paste(nrow(wgs_secondary), "WGS secondary specimens"))
readr::write_tsv(wgs_secondary, wgs_secondary_file)

message(paste(nrow(wgs_primary_plus), "WGS specimens (including non-primary)"))
readr::write_tsv(wgs_primary_plus, wgs_primplus_file)

message(paste(nrow(wgs_primary) + nrow(wxs_panel_primary), "WGS+WXS+Panel primary specimens"))
readr::write_tsv(dplyr::bind_rows(wgs_primary, wxs_panel_primary),
                 wgswxspanel_primary_file)

message(paste(nrow(wgs_secondary) + nrow(wxs_panel_secondary), "WGS+WXS+Panel secondary specimens"))
readr::write_tsv(dplyr::bind_rows(wgs_secondary, wxs_panel_secondary),
                 wgswxspanel_secondary_file)

message(paste(nrow(wgs_primary_plus) + nrow(wxs_panel_primary_plus), "WGS+WXS+Panel specimens (including non-primary)"))
readr::write_tsv(dplyr::bind_rows(wgs_primary_plus, wxs_panel_primary_plus),
                 wgswxspanel_primplus_file)
