# 02-generate-independent-rnaseq.R
# Krutika Gaonkar for D3b
#
# Purpose: Generate tables of independent rna-seq specimens 
# Option descriptions
# -f, --histology_file : File path to where you would like the annotation_rds file to be
#               stored
# -o,--output_directory : Output directory
# example invocation:
# Rscript analyses/independent-samples/02-generate-independent-rnaseq.R \
#   -f data/histologies.tsv \
#   -o analyses/independent-samples/results

# Load the libraries
library(optparse)
library(tidyverse)

# Base directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "independent-samples")

# source sample selection function
source(file.path(analysis_dir, "independent_rna_samples.R"))

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
  ),
  make_option(
    c("-i","--independent_dna_sample_df"),
    type = "character",
    default = NULL,
    help = "path to independent-specimens.wgs* file"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# set output files
out_dir <- opts$output_directory
if (!dir.exists(out_dir)){
  dir.create(out_dir, recursive = TRUE)
}

rnaseq_primplus_file <- file.path(out_dir, 
                                  "independent-specimens.rnaseq.primary-plus.tsv")


# Read histology file
sample_df <- readr::read_tsv(opts$histology_file, 
                             guess_max = 100000,
                             col_types = readr::cols()) # suppress parse message

# Read in dna independent sample list to match to rna samples
# So that independent RNA samples match the DNA samples
independent_dna_sample_df <- read_tsv(opts$independent_dna_sample_df)

# Separating polya and stranded samples since we might have cases where
# we have polya and stranded samples per Kids_First_Participant_ID

# Filter to only samples from tumors, where composition is known to be Solid Tissue
# for all RNA samples
independent_rna_primary_plus <- sample_df %>%
  filter(sample_type == "Tumor", 
         composition == "Solid Tissue"
  ) %>%
  independent_rna_samples(independent_dna_sample_df = 
                            independent_dna_sample_df,
                          histology_df = .,
                          match_type = "independent_dna_plus_only_rna",
                          tumor_description_rna_only = "primary_plus",seed = 2020) %>%
  readr::write_tsv(rnaseq_primplus_file)

