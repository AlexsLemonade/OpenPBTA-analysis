# Authors: Komal S. Rathi and Stephanie J. Spielman
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
              help = "Deconv Output (.RDS)"), 
  make_option("--apply_tumor_purity_threshold",
              action = "store_true",
              default = FALSE,
              help = "This flag turns on filtering by tumor purity. Only stranded data will be processed.")
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

# We read in the stranded data no matter what 
stranded <- readRDS(stranded)


# Filter data if tumor purity threshold is turned on and 
#  only process the relevant stranded libraries (skip polyA)
if (opt$apply_tumor_purity_threshold) {
  
  # tumor purity metadata file which has been filtered to only biospecimens that
  #  survive the cancer-group-level threshold
  tumor_purity_file <- file.path(
    rprojroot::find_root(rprojroot::has_dir(".git")),
    "analyses",
    "tumor-purity-exploration",
    "results",
    "thresholded_rna_stranded_same-extraction.tsv"
  )
  
  # Define vector of IDs to retain
  bs_ids <- readr::read_tsv(tumor_purity_file) %>%
    dplyr::pull(Kids_First_Biospecimen_ID) %>%
    unique()
  
  # Filter the matrix
  stranded <- stranded[, bs_ids]
  
  # Quick check on filtering:
  if (!length(bs_ids) == ncol(stranded)) {
    stop("Error filtering IDs to use in immune deconvolution with tumor purity thresholding on.")
  }
} else {
  
  # If we are not filtering by tumor purity, also read and run polyA
  
  # Read expression data 
  polya <- readRDS(polya)
  
  # Deconvolute polya and add columns to indicate library and method
  result_polya <- deconvolute(gene_expression = as.matrix(polya), method = deconv_method)
  polya_wide <- result_polya %>%
    gather(-cell_type, key = "sample", value = "score") %>%
    mutate(library = "polya", 
           method = deconv_method)
  
}

result_stranded <- deconvolute(gene_expression = as.matrix(stranded), method = deconv_method) %>%
  gather(-cell_type, key = "sample", value = "score") %>%
  # library and method columns
  mutate(library = "stranded",
         method = deconv_method)

# If we are filtering on tumor purity, write out the results
# Otherwise, combine with polyA and then write out results
if (opt$apply_tumor_purity_threshold) {
  write_rds(result_stranded, output_file)
} else {
  combined_results <- result_stranded %>%
    bind_rows(polya_wide) 
  write_rds(combined_results, output_file)
}
