# Author: Komal S. Rathi and Stephanie J. Spielman
# Function:
# Script to perform immune characterization using immunedeconv, uses xCell by default.

# Find the root directory of this repository
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(immunedeconv))

option_list <- list(
  make_option(c("-p", "--polyaexprs"), type = "character",
              help = "PolyA Expression data: HUGO symbol x Sample identifiers (.RDS)"),
  make_option(c("-s", "--strandedexprs"), type = "character",
              help = "Stranded Expression data: HUGO symbol x Sample identifiers (.RDS)"),
  make_option(c("-c", "--clin"), type = "character",
              help = "Clinical file (.TSV)"),
  make_option(c("-m", "--method"), type = "character",
              default = "xcell",
              help = "Deconvolution Method"),
  make_option(c("-o","--outputfile"), type = "character",
              help = "Deconv Output (.RDS)")
)

# Example Run:
# Rscript 01-immune-deconv.R \
# --polyaexprs '../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds' \
# --strandedexprs '../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds' \
# --method 'xcell' \
# --outputfile 'results/xcell_deconv-output.RDS'

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
polya <- opt$polyaexprs
stranded <- opt$strandedexprs
deconv_method <- tolower(opt$method) # immunedeconv requires lower case
output_file <- opt$outputfile

#### Check model parameter - must be in deconvolution_methods (an immmunedeconv variable)
if (!is.null(deconv_method)){
  if (!(deconv_method %in% deconvolution_methods)) {
    stop( paste(c("Specified method not available. Must be one of the following: ", deconvolution_methods), collapse=" ") )
  }
}

# Read expression data 
polya <- readRDS(polya)
stranded <- readRDS(stranded)

# Deconvolute separately and combine results, with columns to indicate method and library
result_polya <- deconvolute(gene_expression = as.matrix(polya), method = deconv_method)
polya_wide <- result_polya %>%
  gather(-cell_type, key = "sample", value = "score") %>%
  mutate(library = "polya")

result_stranded <- deconvolute(gene_expression = as.matrix(stranded), method = deconv_method)

combined_results <- result_stranded %>%
  gather(-cell_type, key = "sample", value = "score") %>%
  mutate(library = "stranded") %>%
  bind_rows(polya_wide) %>%
  mutate(method = deconv_method)

# save output to RDS file
write_rds(combined_results, output_file)
