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
polya.file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.polya.rds")
stranded.file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")
clin.file <- file.path(data_dir, "pbta-histologies.tsv")

# output files
uncorrected.file <- file.path(output_dir, paste0(output_prefix, ".rds"))
corrected.file <- file.path(output_dir, paste0(output_prefix, "-batch-corrected.rds"))

# expression from polya and stranded data
polya <- readRDS(polya.file)
stranded <- readRDS(stranded.file)

# read and subset clinical file to MB samples
clin <- read.delim(clin.file, stringsAsFactors = F)
clin.mb  <- clin %>%
  filter(experimental_strategy == "RNA-Seq",
         integrated_diagnosis == "Medulloblastoma") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, RNA_library)

# function to filter only MB samples
filter.mat <- function(expr.input, clin.mb) {
  
  # get data
  expr.input <- get(expr.input)
  
  # subset expression to MB samples
  mb.samples <- intersect(clin.mb$Kids_First_Biospecimen_ID, colnames(expr.input))
  expr.input <- expr.input %>%
    dplyr::select(mb.samples)
  
  return(expr.input)
}

# filter matrices to MB only 
expr.input <- c("polya", "stranded") 
expr.input.mb <- lapply(expr.input, FUN = function(x) filter.mat(expr.input = x, clin = clin.mb))

# combine both matrices using common genes 
expr.input.mb <- expr.input.mb[[1]] %>%
  rownames_to_column('gene') %>%
  inner_join(expr.input.mb[[2]] %>%
               rownames_to_column('gene'), by = 'gene') %>%
  column_to_rownames('gene')

# save uncorrected matrix
uncorrected.df <- log2(expr.input.mb + 1)
write_rds(uncorrected.df, uncorrected.file)

# batch correct if batch_col not null
if(!is.null(batch_col)){
  print("Batch correct input matrices...")
  
  # match clinical rows and expression cols
  expr.input.mb  <- as.matrix(expr.input.mb[,clin.mb$Kids_First_Biospecimen_ID])
  
  # batch correct using batch_col
  corrected.mat <- ComBat(dat = log2(expr.input.mb + 1), batch = clin.mb[, batch_col])
  
  # save corrected matrix
  write_rds(corrected.mat, corrected.file)
}

