# Author: Komal S. Rathi
# Date: 11/11/2019
# Function: 
# Script to perform immune characterization using xCell etc.

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
  make_option(c("-b", "--cibersortbin"), type = "character",
              help = "Path to Cibersort binary (CIBERSORT.R)"),
  make_option(c("-m", "--cibersortmat"), type = "character",
              help = "Path to Cibersort matrix (LM22.txt)"),
  make_option(c("-o","--outputfile"), type = "character",
              help = "Deconv Output (.RData)")
)

# Example Run:
# Rscript analyses/immune-deconv/01-immune-deconv.R \
# -p 'data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds' \
# -s 'data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds' \
# -c 'data/pbta-histologies.tsv' \
# -b 'analyses/immune-deconv/CIBERSORT.R' \
# -m 'analyses/immune-deconv/LM22.txt' \
# -o 'analyses/immune-deconv/deconv-output.RData'

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
polya <- opt$polyaexprs
stranded <- opt$strandedexprs
clin.file <- opt$clin
output.file <- opt$outputfile

# set path to cibersort binary and matrix
cibersort_bin <- opt$cibersortbin 
set_cibersort_binary(cibersort_bin)
cibersort_mat <- opt$cibersortmat 
set_cibersort_mat(cibersort_mat)

# merge expression from polya and stranded data on common genes
polya <- readRDS(polya)
stranded <- readRDS(stranded)
common.genes <- intersect(rownames(polya), rownames(stranded))
polya <- polya[common.genes,]
stranded <- stranded[common.genes,]
expr.input <- cbind(polya, stranded)

# read clinical data
clin <- read.delim(clin.file, stringsAsFactors = F)
clin  <- clin %>% 
  filter(Kids_First_Biospecimen_ID %in% colnames(expr.input)) %>%
  select(Kids_First_Biospecimen_ID, broad_histology)

# function to run immunedeconv
deconv <- function(expr.input, method){
  
  # deconvolute using specified method
  res <- deconvolute(gene_expression = as.matrix(expr.input), method = method)
  res$method <- method # assign method name 

  # merge output with clinical data
  res <- res %>%
    gather(sample, fraction, -c(cell_type, method)) %>%
    as.data.frame()
  res <- merge(res, clin, by.x = 'sample', by.y = 'Kids_First_Biospecimen_ID')

  return(res)
}
  
# Deconvolute using xCell and Cibersort (absolute mode)
# these two methods have the max number of immune cell types
deconv.methods <- immunedeconv::deconvolution_methods
res <- sapply(deconv.methods[grep("xcell|cibersort_abs", deconv.methods)], function(x) deconv(expr.input, method = x))
xcell <- as.data.frame(res[,'xCell'])
cibersort_abs <- as.data.frame(res[,'CIBERSORT (abs.)'])

# save all deconvolution methods' output to one RData object 
print("Writing output to file..")
save(cibersort_abs, xcell, file = output.file)
print("Done!")
