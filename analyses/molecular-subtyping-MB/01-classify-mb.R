# Author: Komal S. Rathi
# Date: 07/22/2020
# Function:
# Script to perform MB molecular subtyping

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
  make_option(c("--corrected_mat"), type = "character",
              help = "PolyA Expression data: HUGO symbol x Sample identifiers (.rds)"),
  make_option(c("--uncorrected_mat"), type = "character",
              help = "Stranded Expression data: HUGO symbol x Sample identifiers (.rds)"),
  make_option(c("--output_prefix"), type = "character",
              help = "Output file prefix")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
corrected_mat <- opt$corrected_mat
uncorrected_mat <- opt$uncorrected_mat
method <- opt$method
output_prefix <- opt$output_prefix

# expression from polya and stranded data
corrected_mat <- readRDS(corrected_mat)
uncorrected_mat <- readRDS(uncorrected_mat)

# combination of dataset and methods
methods <- c('MM2S', 'medullo-classifier')
mats <- c('corrected_mat', 'uncorrected_mat')
combo <- expand.grid(mats, methods, stringsAsFactors = F) 

# classify mb samples
print("Classify medulloblastoma subtypes...")
mb.classify <- apply(combo, 1, FUN = function(x) classify.mb(expr.input = x[1], method = x[2]))
names(mb.classify) <- apply(combo, 1, FUN = function(x) paste0(x[1], '_', x[2]))

# save output to rds object
print("Writing output to file..")
outputfile <- file.path(output_dir, paste0(output_prefix, ".rds"))
saveRDS(mb.classify, file = outputfile)
print("Done!")
