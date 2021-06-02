# J. Taroni for ALSF CCDL 2020
# This script accepts the OpenPBTA histologies file + a SEG file and will
# produce an array list file that contains WGS samples that have the
# user-specified filter_value in the user-specified filter_column.
# We include the SEG file so we can take the intersection of the IDs that match
# the criteria
#
# Example usage - filtering to HGAT samples
# (assumes you are in the analyses/run-gistic directory)
#
#   Rscript --vanilla scripts/subset-seg-file.R \
#     --segfile ../../data/pbta-cnv-consensus.seg.gz \
#     --metadata ../../data/pbta-histologies.tsv \
#     --filter_column short_histology \
#     --filter_value HGAT \
#     --output_file array_list_files/hgat-array-file.txt
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

#### Generate array list file --------------------------------------------------

relevant_biospecimen_ids <- histologies_df %>%
  # only WGS samples have copy number data
  filter(experimental_strategy == "WGS",
         # ID rows that contain the filter_value in filter_column
         !!sym(opt$filter_column) == opt$filter_value) %>%
  pull(Kids_First_Biospecimen_ID)

# Get the intersection of the IDs in the SEG file with the IDs that meet the
# criteria -- some samples may have been filtered out during the consensus call
# process for example
relevant_biospecimen_ids <- intersect(seg_df$ID,
                                      relevant_biospecimen_ids)

# write to file - it is a text file with a single column with an optional
# header "array"
write_delim(data.frame(array = relevant_biospecimen_ids),
            path = opt$output_file)
