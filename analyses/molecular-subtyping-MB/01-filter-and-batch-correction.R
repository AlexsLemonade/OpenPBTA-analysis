# Author: Komal S. Rathi
# Date: 08/19/2020
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
clin_file <- file.path(data_dir, "pbta-histologies-base.tsv")
polya.file <- file.path(root_dir, "data", "pbta-gene-expression-rsem-fpkm-collapsed.polya.rds")
stranded.file <- file.path(root_dir, "data","pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")
terms.file <- file.path("input", "mb_subtyping_path_dx_strings.json")

# Read in the JSON file that contains the strings we'll use to include or
# exclude samples for subtyping - see 00-mb-select-pathology-dx
path_dx_list <- jsonlite::fromJSON(terms.file)

# output files
uncorrected.file <- file.path(output_dir, paste0(output_prefix, ".rds"))
corrected.file <- file.path(output_dir, paste0(output_prefix, "-batch-corrected.rds"))

# expression from polya and stranded data
polya <- readRDS(polya.file)
stranded <- readRDS(stranded.file)

# read and subset clinical file to MB samples
clin <- read.delim(clin_file, stringsAsFactors = F)

# Filter to medulloblastoma samples only based on criteria determined in 00-mb-select-pathology-dx
clin.mb  <- clin %>%
  # Inclusion on the basis of strings in pathology_diagnosis and
  # if 'Other', then based on pathology_free_text_diagnosis
  filter(pathology_diagnosis %in% path_dx_list$exact_path_dx |
         (pathology_diagnosis == "Other" & str_detect(str_to_lower(pathology_free_text_diagnosis), paste0(path_dx_list$include_free_text, collapse = "|"))))

# Write to TSV for use later
readr::write_tsv(clin.mb, file.path("input", "subset-mb-clinical.tsv"))

# Keep to only RNA-seq for these next step
clin.mb.rnaseq <- clin.mb %>%
  dplyr::filter(experimental_strategy == "RNA-Seq")

# Make this a expression so we can match any of these
biospecimen_ids <- paste0(clin.mb.rnaseq$Kids_First_Biospecimen_ID, collapse = "|")

# Select only the mb biospecimens in polya
polya.mb <- polya %>%
  dplyr::select(dplyr::matches(biospecimen_ids)) %>%
  rownames_to_column('gene')

# Select only the mb biospecimens in stranded
stranded.mb <-  stranded %>%
  dplyr::select(dplyr::matches(biospecimen_ids)) %>%
  rownames_to_column('gene')

# inner_join by gene, turn into a matrix
all.exprs.mb <- polya.mb %>%
  inner_join(stranded.mb, by = "gene") %>%
  column_to_rownames('gene') %>%
  as.matrix()

# save uncorrected matrix
uncorrected.mat <- log2(all.exprs.mb + 1)
write_rds(uncorrected.mat, uncorrected.file)

# batch correct if batch_col not null
if(!is.null(batch_col)){
  print("Batch correct input matrices...")

  # match clinical rows and expression cols
  expr.input.mb  <- as.matrix(all.exprs.mb[, clin.mb.rnaseq$Kids_First_Biospecimen_ID])

  # batch correct using batch_col
  corrected.mat <- ComBat(dat = log2(expr.input.mb + 1), batch = clin.mb.rnaseq[, batch_col])

  # save corrected matrix
  write_rds(corrected.mat, corrected.file)
}
