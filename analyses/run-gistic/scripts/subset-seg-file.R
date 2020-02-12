# J. Taroni for ALSF CCDL 2020
# This script accepts the OpenPBTA histologies file + a SEG file and will filter
# the SEG file to the WGS samples that have the user-specified filter_value in
# the user-specified filter_column
#
# Example usage - filtering to HGAT samples
# (assumes you are in the analyses/run-gistic directory)
#
#   Rscript --vanilla scripts/subset-seg-file.R \
#     --segfile ../../data/pbta-cnv-consensus.seg.gz \
#     --metadata ../../data/pbta-histologies.tsv \
#     --filter_column short_histology \
#     --filter_value HGAT \
#     --output_file seg_files/hgat-cnv-consensus.seg.gz
#

library(tidyverse)

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("--segfile"),
    type = "character",
    default = NULL,
    help = "Full path SEG file to subset",
  ),
  optparse::make_option(
    c("--metadata"),
    type = "character",
    default = NULL,
    help = "Full path to metadata TSV file"
  ),
  optparse::make_option(
    c("--filter_column"),
    type = "character",
    default = "short_histology",
    help = "Name of the column that we will filter on"
  ),
  optparse::make_option(
    c("--filter_value"),
    type = "character",
    default = NULL,
    help = "Value to match in the filter_column to qualify for inclusion"
  ),
  optparse::make_option(
    c("--output_file"),
    type = "character",
    default = NULL,
    help = "A character vector that will be used to name output files"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#### Read in files -------------------------------------------------------------

seg_df <- read_tsv(opt$segfile)
histologies_df <- read_tsv(opt$metadata)

#### Filter SEG file -----------------------------------------------------------

relevant_biospecimen_ids <- histologies_df %>%
  # only WGS samples have copy number data
  filter(experimental_strategy == "WGS",
         # ID rows that contain the filter_value in filter_column
         !!sym(opt$filter_column) == opt$filter_value) %>%
  pull(Kids_First_Biospecimen_ID)

seg_df %>%
  filter(ID %in% relevant_biospecimen_ids) %>%
  write_tsv(opt$output_file)
