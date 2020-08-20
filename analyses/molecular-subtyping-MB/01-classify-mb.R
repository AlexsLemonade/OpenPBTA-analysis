# Author: Komal S. Rathi
# Date: 07/22/2020
# Function:
# Script to perform MB molecular subtyping

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(medulloPackage))
suppressPackageStartupMessages(library(MM2S))
suppressPackageStartupMessages(library(org.Hs.eg.db))

# source classification function
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "molecular-subtyping-MB", "util", "classify-mb.R"))

option_list <- list(
  make_option(c("--polyaexprs"), type = "character",
              help = "PolyA Expression data: HUGO symbol x Sample identifiers (.rds)"),
  make_option(c("--strandedexprs"), type = "character",
              help = "Stranded Expression data: HUGO symbol x Sample identifiers (.rds)"),
  make_option(c("--batch_col"), type = "character", 
              default = NULL,
              help = "Combine and batch correct input matrices using which column?"),
  make_option(c("--clin"), type = "character",
              help = "Clinical file (.tsv)"),
  make_option(c("--hist_column"), type = "character",
              help = "Histology column to subset clinical file to Medulloblastoma samples"),
  make_option(c("--method"), type = "character", 
              help = "Method to use for classification: MM2S or medullo-classifier"),
  make_option(c("--outputfile"), type = "character",
              help = "Output file path (.rds)")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
polya <- opt$polyaexprs
stranded <- opt$strandedexprs
batch_col <- opt$batch_col
clin.file <- opt$clin
hist_column <- opt$hist_column
method <- opt$method
outputfile <- opt$outputfile

# check input method 
if(!method %in% c("MM2S", "medullo-classifier")){
  stop("Check method name")
}

# expression from polya and stranded data
polya <- readRDS(polya)
stranded <- readRDS(stranded)

# read and subset clinical file to MB samples
clin <- read.delim(clin.file, stringsAsFactors = F)
clin.mb  <- clin %>%
  filter(experimental_strategy == "RNA-Seq",
         !!as.name(hist_column) == "Medulloblastoma") %>%
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

# if batch_col is set, do the following:
if(!is.null(batch_col)){
  print("Batch correct input matrices...")
  
  # match clinical rows and expression cols
  expr.input.mb  <- as.matrix(expr.input.mb[,clin.mb$Kids_First_Biospecimen_ID])
  
  # batch correct using batch_col
  expr.input.mb <- ComBat(dat = log2(expr.input.mb + 1), batch = clin.mb[, batch_col])
} else {
  # log transform
  expr.input.mb <- log2(expr.input.mb + 1)
}

# classify mb samples
print("Classify medulloblastoma subtypes...")
mb.classify <- classify.mb(expr.mb = expr.input.mb, clin.mb = clin.mb, method = method)

# save output to RData object
print("Writing output to file..")
saveRDS(mb.classify, file = outputfile)
print("Done!")
