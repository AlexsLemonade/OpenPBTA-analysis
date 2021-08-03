# Author: Komal S. Rathi
# Function: Script to filter MB samples and/or batch correct

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sva))

option_list <- list(
  make_option(c("--batch_col"), type = "character",
              default = NULL,
              help = "Combine and batch correct input matrices using which column?"),
  make_option(c("--output_dir"), type = "character",
              help = "Output directory"),
  make_option(c("--output_prefix"), type = "character",
              help = "Output file prefix")
)

# parse options
opt <- parse_args(OptionParser(option_list = option_list))
batch_col <- opt$batch_col
output_prefix <- opt$output_prefix
output_dir <- opt$output_dir

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# input files
clin_file <- file.path(data_dir, "histologies.tsv")
exprs_file <- file.path(data_dir, "gene-expression-rsem-tpm-collapsed.rds")
terms_file <- file.path("input", "mb_subtyping_path_dx_strings.json")

# Read in the JSON file that contains the strings we'll use to include or
# exclude samples for subtyping - see 00-mb-select-pathology-dx
path_dx_list <- jsonlite::fromJSON(terms_file)

# output file
uncorrected_file <- file.path(output_dir, paste0(output_prefix, ".rds"))

# collapsed expression matrix
exprs_mat <- readRDS(exprs_file)

# read and subset clinical file to MB samples
clin <- read.delim(clin_file, stringsAsFactors = F)

# Filter to medulloblastoma samples only based on criteria determined in 00-mb-select-pathology-dx
clin_mb  <- clin %>%
  # Inclusion on the basis of strings in pathology_diagnosis and
  # if 'Other', then based on pathology_free_text_diagnosis
  filter(pathology_diagnosis %in% path_dx_list$exact_path_dx |
           (pathology_diagnosis == "Other" & str_detect(str_to_lower(pathology_free_text_diagnosis), paste0(path_dx_list$include_free_text, collapse = "|"))))

# Write to TSV for use later
readr::write_tsv(clin_mb, file.path("input", "subset-mb-clinical.tsv"))

# Keep to only RNA-seq for these next step
clin_mb_rnaseq <- clin_mb %>%
  dplyr::filter(experimental_strategy == "RNA-Seq")

# Select only the mb biospecimens in expression mat
exprs_mat_mb <- exprs_mat %>%
  dplyr::select(clin_mb_rnaseq$Kids_First_Biospecimen_ID)

# save uncorrected matrix
uncorrected_mat <- log2(exprs_mat_mb + 1)
write_rds(uncorrected_mat, uncorrected_file)

# batch correct if batch_col not null
if(!is.null(batch_col)){

  print("Batch correct input matrices...")
  
  # batch correct using batch_col
  corrected_mat <- ComBat(dat = log2(exprs_mat_mb + 1), 
                          batch = clin_mb_rnaseq[, batch_col])
  
  # save corrected matrix
  corrected_file <- file.path(output_dir, paste0(output_prefix, "-batch-corrected.rds"))
  write_rds(corrected_mat, corrected_file)
}
