# Author: Komal S. Rathi
# Date: 11/11/2019
# Function:
# Script to perform immune characterization using xCell etc.

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
              help = "Deconvolution Method"),
  make_option(c("-b", "--cibersortbin"), type = "character",
              help = "Path to Cibersort binary (CIBERSORT.R)"),
  make_option(c("-g", "--cibersortgenemat"), type = "character",
              help = "Path to Cibersort signature matrix (LM22.txt)"),
  make_option(c("-o","--outputfile"), type = "character",
              help = "Deconv Output (.RData)")
)

# Example Run:
# Rscript analyses/immune-deconv/01-immune-deconv.R \
# -p 'data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds' \
# -s 'data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds' \
# -c 'data/pbta-histologies.tsv' \
# -m $DECONV_METHOD \
# -b $CIBERSORT_BIN \
# -g $CIBERSORT_MAT \
# -o 'analyses/immune-deconv/results/deconv-output.RData'

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
polya <- opt$polyaexprs
stranded <- opt$strandedexprs
clin.file <- opt$clin
deconv.method <- opt$method
cibersort_bin <- opt$cibersortbin
cibersort_mat <- opt$cibersortgenemat
output.file <- opt$outputfile

#### Check model parameter - must be in deconvolution_methods (immunedeconv accepted options)
if (!is.null(deconv.method)){
  if (!(deconv.method %in% deconvolution_methods)) {
    stop( paste(c("Specified method not available. Must be one of the following: ", deconvolution_methods), collapse=" ") )
  }
}


# if cibersort_bin and cibersort_mat are defined
# then, set path to cibersort binary and matrix
if(!(is.null(cibersort_bin)) & !(is.null(cibersort_mat))){
  set_cibersort_binary(cibersort_bin)
  set_cibersort_mat(cibersort_mat)
}

# merge expression from polya and stranded data on common genes
polya <- readRDS(polya)
stranded <- readRDS(stranded)

# read clinical data
clin <- readr::read_tsv(clin.file, guess_max = 10000)

# function to run immunedeconv
deconv <- function(expr.input, method) {

  # get data
  expr.input <- get(expr.input)

  # Import standard color palettes for project
  histology_label_mapping <- readr::read_tsv(
    file.path(root_dir, "figures", "palettes", "histology_label_color_table.tsv")) %>%
    # Select just the columns we will need for plotting
    dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes) %>%
    # Reorder display_group based on display_order
    dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))

  # subset clinical
  clin.sub  <- clin %>%
    filter(Kids_First_Biospecimen_ID %in% colnames(expr.input)) %>%
    dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID") %>%
    dplyr::select(Kids_First_Biospecimen_ID, broad_histology, display_group, molecular_subtype)

  # deconvolute using specified method
  res <- deconvolute(gene_expression = as.matrix(expr.input), method = method)
  res$method <- names(grep(method, deconvolution_methods, value = TRUE)) # assign method name

  # merge output with clinical data
  res <- res %>%
    gather(sample, fraction, -c(cell_type, method)) %>%
    as.data.frame() %>%
    inner_join(clin.sub, by = c("sample" = "Kids_First_Biospecimen_ID"))

  return(res)
}

# Deconvolute using xCell and the second chosen method
deconv.method <- unique(c("xcell", deconv.method))
expr.input <- c("polya", "stranded") # datasets
combo <- expand.grid(expr.input, deconv.method, stringsAsFactors = F) # combination of dataset and methods
deconv.res <- apply(combo, 1, FUN = function(x) deconv(expr.input = x[1], method = x[2]))
deconv.res <- do.call(rbind.data.frame, deconv.res)

# save output to RData object
print("Writing output to file..")
save(deconv.res, file = output.file)
print("Done!")
