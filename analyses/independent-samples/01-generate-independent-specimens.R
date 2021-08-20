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
source(file.path(analysis_dir, "util", "independent-samples.R"))

set.seed(2020)

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

wgs_primary_each_file <- file.path(out_dir, 
                              "independent-specimens.wgs.primary.eachcohort.tsv")
wgs_relapse_each_file <- file.path(out_dir, 
                              "independent-specimens.wgs.relapse.eachcohort.tsv")
wgs_primplus_each_file <- file.path(out_dir, 
                              "independent-specimens.wgs.primary-plus.eachcohort.tsv")
wgswxspanel_primary_each_file <- file.path(out_dir, 
                              "independent-specimens.wgswxspanel.primary.eachcohort.tsv")
wgswxspanel_relapse_each_file <- file.path(out_dir, 
                              "independent-specimens.wgswxspanel.relapse.eachcohort.tsv")
wgswxspanel_primplus_each_file <- file.path(out_dir, 
                              "independent-specimens.wgswxspanel.primary-plus.eachcohort.tsv")

wgs_primary_all_file <- file.path(out_dir, 
                                   "independent-specimens.wgs.primary.tsv")
wgs_relapse_all_file <- file.path(out_dir, 
                                   "independent-specimens.wgs.relapse.tsv")
wgs_primplus_all_file <- file.path(out_dir, 
                                    "independent-specimens.wgs.primary-plus.tsv")
wgswxspanel_primary_all_file <- file.path(out_dir, 
                                           "independent-specimens.wgswxspanel.primary.tsv")
wgswxspanel_relapse_all_file <- file.path(out_dir, 
                                           "independent-specimens.wgswxspanel.relapse.tsv")
wgswxspanel_primplus_all_file <- file.path(out_dir, 
                                            "independent-specimens.wgswxspanel.primary-plus.tsv")

# Read histology file
histology_df <- readr::read_tsv(opts$histology_file, 
                             guess_max = 100000,
                             col_types = readr::cols()) # suppress parse message

# randomize rows of histology file to avoid selection bias
set.seed(100)
histology_df <- histology_df[sample(nrow(histology_df)), ]

# Filter to only samples from tumors, where composition is known to be Solid Tissue or Bone Marrow
# Note that there are some samples with unknown composition, but these will be ignored for now.
tumor_samples <- histology_df %>%
  dplyr::filter(sample_type == "Tumor", 
                composition == "Solid Tissue" | composition == "Bone Marrow", 
                experimental_strategy %in% c("WGS", "WXS", "Targeted Sequencing"))

# Generate WGS independent samples
wgs_samples <- tumor_samples %>%
  dplyr::filter(experimental_strategy == "WGS")

wgs_primary_each <- independent_samples(wgs_samples, tumor_types = "primary", independent_level = "each-cohort", seed = 2020)
wgs_relapse_each <- independent_samples(wgs_samples, tumor_types = "relapse", independent_level = "each-cohort", seed = 2020)
wgs_primary_plus_each <- independent_samples(wgs_samples, tumor_types = "prefer_primary", independent_level = "each-cohort", seed = 2020)

wgs_primary_all <- independent_samples(wgs_samples, tumor_types = "primary", independent_level = "all-cohorts", seed = 2020)
wgs_relapse_all <- independent_samples(wgs_samples, tumor_types = "relapse", independent_level = "all-cohorts", seed = 2020)
wgs_primary_plus_all <- independent_samples(wgs_samples, tumor_types = "prefer_primary", independent_level = "all-cohorts", seed = 2020)

# Generate lists for WXS and Panel samples 
# WGS is generally preferred, so we will only include those where WGS is not available
wxs_panel_samples <-  tumor_samples %>% 
  dplyr::filter(!(Kids_First_Participant_ID %in% 
                  wgs_samples$Kids_First_Participant_ID))

wxs_panel_primary_each <- independent_samples(wxs_panel_samples, tumor_types = "primary", independent_level = "each-cohort", seed = 2020)
wxs_panel_relapse_each <- independent_samples(wxs_panel_samples, tumor_types = "relapse", independent_level = "each-cohort", seed = 2020)
wxs_panel_primary_plus_each <- independent_samples(wxs_panel_samples, tumor_types = "prefer_primary", independent_level = "each-cohort", seed = 2020)

wxs_panel_primary_all <- independent_samples(wxs_panel_samples, tumor_types = "primary", independent_level = "all-cohorts", seed = 2020)
wxs_panel_relapse_all <- independent_samples(wxs_panel_samples, tumor_types = "relapse", independent_level = "all-cohorts", seed = 2020)
wxs_panel_primary_plus_all <- independent_samples(wxs_panel_samples, tumor_types = "prefer_primary", independent_level = "all-cohorts", seed = 2020)

# write files for independent specimens considering cohort difference - for WGS specimens only 
message(paste(nrow(wgs_primary_each), "WGS primary specimens for each cohort"))
wgs_primary_each %>% dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_tsv(wgs_primary_each_file)

message(paste(nrow(wgs_relapse_each), "WGS relapse specimens for each cohort"))
wgs_relapse_each %>% dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgs_relapse_each_file)

message(paste(nrow(wgs_primary_plus_each), "WGS specimens (including non-primary) for each cohort"))
wgs_primary_plus_each %>% dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgs_primplus_each_file)

# write files for independent specimens not considering cohort difference
message(paste(nrow(wgs_primary_all), "WGS primary specimens for all cohorts"))
wgs_primary_all %>% dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgs_primary_all_file)

message(paste(nrow(wgs_relapse_all), "WGS relapse specimens for all cohorts"))
wgs_relapse_all %>% dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgs_relapse_all_file)

message(paste(nrow(wgs_primary_plus_all), "WGS specimens (including non-primary) for all cohorts"))
wgs_primary_plus_all %>% dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgs_primplus_all_file)

# write files for independent specimens considering cohort difference - for WGS+WXS+Panel specimens 
message(paste(nrow(wgs_primary_each) + nrow(wxs_panel_primary_each), "WGS+WXS+Panel primary specimens for each cohort"))
dplyr::bind_rows(wgs_primary_each, wxs_panel_primary_each) %>% 
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgswxspanel_primary_each_file)

message(paste(nrow(wgs_relapse_each) + nrow(wxs_panel_relapse_each), "WGS+WXS+Panel relapse specimens for each cohort"))
dplyr::bind_rows(wgs_relapse_each, wxs_panel_relapse_each) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgswxspanel_relapse_each_file)

message(paste(nrow(wgs_primary_plus_each) + nrow(wxs_panel_primary_plus_each), "WGS+WXS+Panel specimens (including non-primary) for each cohort"))
dplyr::bind_rows(wgs_primary_plus_each, wxs_panel_primary_plus_each) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgswxspanel_primplus_each_file)

# write files for independent specimens not considering cohort difference - for WGS+WXS+Panel specimens 

message(paste(nrow(wgs_primary_all) + nrow(wxs_panel_primary_all), "WGS+WXS+Panel primary specimens for all cohort"))
dplyr::bind_rows(wgs_primary_all, wxs_panel_primary_all) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgswxspanel_primary_all_file)

message(paste(nrow(wgs_relapse_all) + nrow(wxs_panel_relapse_all), "WGS+WXS+Panel relapse specimens for all cohort"))
dplyr::bind_rows(wgs_relapse_all, wxs_panel_relapse_all) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgswxspanel_relapse_all_file)

message(paste(nrow(wgs_primary_plus_all) + nrow(wxs_panel_primary_plus_all), "WGS+WXS+Panel specimens (including non-primary) for all cohort"))
dplyr::bind_rows(wgs_primary_plus_all, wxs_panel_primary_plus_all) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(wgswxspanel_primplus_all_file)


