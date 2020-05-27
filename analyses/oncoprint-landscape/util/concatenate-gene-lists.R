# This script concatenates given gene lists for the purpose of plotting.
#
# Chante Bethell for CCDL 2020
#
# EXAMPLE USAGE:
#
# Rscript --vanilla concatenate-gene-lists.R \
#  --goi_file_1 ../interaction-plots/scratch/all_primary_samples_maf.tsv \
#  --goi_file_2 ../../scratch/all_primary_samples_cnv.tsv

#### Set Up --------------------------------------------------------------------

# Load tidyverse
library(tidyverse)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define file path to output directory
results_dir <-
  file.path(root_dir, "analyses", "oncoprint-landscape", "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-a", "--goi_file_1"),
    type = "character",
    default = NULL,
    help = "file path to first file with genes of interest",
  ),
  optparse::make_option(
    c("-b", "--goi_file_2"),
    type = "character",
    default = NULL,
    help = "file path to second file with genes of interest"
  ),
  optparse::make_option(
    c("-n", "--filename"),
    type = "character",
    default = NULL,
    help = "string for the output filename"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Error handling related to specifying a second goi_list
if (is.null(opt$goi_file_2)) {
  stop("You must specify the file path to second genes of interest file")
}

#### Read in data --------------------------------------------------------------

# Read in genes of interest file 1 and isolate the genes
goi_df_1 <-
  read_tsv(file.path(opt$goi_file_1)) %>%
  select(gene)

# Read in genes of interest file 2 and isolate the genes
goi_df_2 <-
  read_tsv(file.path(opt$goi_file_2)) %>%
  select(gene)

#### Concatenate gene lists ---------------------------------------------------

final_gene_list <- goi_df_1 %>%
  bind_rows(goi_df_2) %>%
  # Get rid of duplicate genes that result from an overlap between the gene lists
  # and removal of extraneous variables
  distinct() %>%
  # Write to file
  write_tsv(file.path(results_dir, opt$filename))
