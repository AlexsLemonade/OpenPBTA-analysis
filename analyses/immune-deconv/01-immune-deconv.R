# Author: Komal S. Rathi
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
              help = "Deconv Output (.RData)")
)

# Example Run:
# Rscript 01-immune-deconv.R \
# --polyaexprs '../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds' \
# --strandedexprs '../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds' \
# --clin '../../data/pbta-histologies.tsv' \
# --method 'xcell' \
# --outputfile 'results/deconv-output.RData'

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
polya <- opt$polyaexprs
stranded <- opt$strandedexprs
clin_file <- opt$clin
deconv_method <- opt$method
output_file <- opt$outputfile

#### Check model parameter - must be in deconvolution_methods (immunedeconv accepted options)
if (!is.null(deconv_method)){
  if (!(deconv_method %in% deconvolution_methods)) {
    stop( paste(c("Specified method not available. Must be one of the following: ", deconvolution_methods), collapse=" ") )
  }
}

# merge expression from polya and stranded data on common genes
polya <- readRDS(polya)
stranded <- readRDS(stranded)

# read clinical data
clin_file <- readr::read_tsv(clin_file, guess_max = 10000)

# function to run immunedeconv
deconv <- function(expr_input, clin_file, method) {

  # get data
  expr_input <- get(expr_input)

  # Import standard color palettes for project
  histology_label_mapping <- readr::read_tsv(
    file.path(root_dir, "figures", "palettes", "histology_label_color_table.tsv")) %>%
    # Select just the columns we will need for plotting
    dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes) %>%
    # Reorder display_group based on display_order
    dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))

  # subset clinical
  clin_file_sub  <- clin_file %>%
    filter(Kids_First_Biospecimen_ID %in% colnames(expr_input)) %>%
    dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID") %>%
    dplyr::select(Kids_First_Biospecimen_ID, broad_histology, display_group, molecular_subtype)

  # deconvolute using specified method
  res <- deconvolute(gene_expression = as.matrix(expr_input), method = method)
  res$method <- names(grep(method, deconvolution_methods, value = TRUE)) # assign method name

  # merge output with clinical data
  res <- res %>%
    gather(sample, fraction, -c(cell_type, method)) %>%
    as.data.frame() %>%
    inner_join(clin_file_sub, by = c("sample" = "Kids_First_Biospecimen_ID"))

  return(res)
}

# Deconvolute using xCell for poly-a and stranded datasets
expr_input <- c("polya", "stranded")
combo <- expand.grid(expr_input, deconv_method, stringsAsFactors = F) 
deconv_output <- apply(combo, 1, FUN = function(x) deconv(expr_input = x[1], clin_file = clin_file, method = x[2]))
deconv_output <- do.call(rbind.data.frame, deconv_output)

# save output to RData object
print("Writing output to file..")
save(deconv_output, file = output_file)
print("Done!")
