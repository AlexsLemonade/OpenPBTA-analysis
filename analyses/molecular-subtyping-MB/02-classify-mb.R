# Author: Komal S. Rathi
# Function: Script to perform MB molecular subtyping

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(medulloPackage))
suppressPackageStartupMessages(library(MM2S))
suppressPackageStartupMessages(library(org.Hs.eg.db))

# source classification function
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "molecular-subtyping-MB", "util", "classify-mb.R"))

# create results directory
output_dir <- file.path(root_dir, "analyses", "molecular-subtyping-MB", "results") 
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

option_list <- list(
  make_option(c("--exprs_mat"), type = "character",
              help = "Expression data: HUGO symbol x Sample identifiers (.rds)"),
  make_option(c("--data_type"), type = "character",
              help = "Expression data type: uncorrected or batch-corrected"),
  make_option(c("--output_prefix"), type = "character",
              help = "Output file prefix")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))

exprs_mat <- opt$exprs_mat
data_type <- opt$data_type
output_prefix <- opt$output_prefix

# expression data
exprs_mat <- readRDS(exprs_mat)

# subtyping methods
methods <- c('MM2S', 'medulloPackage')

# classify mb samples
print("Classify medulloblastoma subtypes...")
classify_mb_output <- lapply(methods, FUN = function(x) classify_mb(exprs_mat = exprs_mat, data_type = data_type, method = x))
names(classify_mb_output) <- methods

# save output to rds object
print("Writing output to file..")
outputfile <- file.path(output_dir, paste0(output_prefix, ".rds"))
saveRDS(classify_mb_output, file = outputfile)
print("Done!")
