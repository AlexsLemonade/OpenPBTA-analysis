# Author: Komal S. Rathi
# Date: 07/22/2020
# Function:
# Script to perform MB molecular subtyping

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(medulloPackage))

option_list <- list(
  make_option(c("--polyaexprs"), type = "character",
              help = "PolyA Expression data: HUGO symbol x Sample identifiers (.rds)"),
  make_option(c("--strandedexprs"), type = "character",
              help = "Stranded Expression data: HUGO symbol x Sample identifiers (.rds)"),
  make_option(c("--batch_col"), type = "character", 
              default = NA,
              help = "Combine and batch correct input matrices using which column?"),
  make_option(c("--clin"), type = "character",
              help = "Clinical file (.tsv)"),
  make_option(c("--hist_column"), type = "character",
              help = "Histology column to subset clinical file to Medulloblastoma samples"),
  make_option(c("--outputfile"), type = "character",
              help = "Subtyping Output (.rds)")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
polya <- opt$polyaexprs
stranded <- opt$strandedexprs
batch_col <- opt$batch_col
clin.file <- opt$clin
output.file <- opt$outputfile
hist_column <- opt$hist_column

# expression from polya and stranded data
polya <- readRDS(polya)
stranded <- readRDS(stranded)

# read and subset clinical file to MB samples
clin <- read.delim(clin.file, stringsAsFactors = F)
clin.mb  <- clin %>%
  filter(experimental_strategy == "RNA-Seq",
         !!as.name(hist_column) == "Medulloblastoma") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, RNA_library)

# filter to MB samples
filter.mat <- function(expr.input, clin.mb) {

  # get data
  expr.input <- get(expr.input)
  
  # subset expression to MB samples
  mb.samples <- intersect(clin.mb$Kids_First_Biospecimen_ID, colnames(expr.input))
  expr.input <- expr.input %>%
    dplyr::select(all_of(mb.samples))
  
  return(expr.input)
}

# function to run molecular subtyping
classify.mat <- function(expr.mb, clin.mb){
  res <- medulloPackage::classify(expr.mb)
  res <- clin.mb %>%
    inner_join(res, by = c("Kids_First_Biospecimen_ID" = "sample"))
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
if(!is.na(batch_col)){
  print("Batch correct input matrices...")
  
  # match clinical rows and expression cols
  clin.mb <- clin.mb %>%
    mutate(tmp = Kids_First_Biospecimen_ID) %>%
    column_to_rownames('tmp')
  expr.input.mb  <- expr.input.mb[,rownames(clin.mb)]
  
  # batch correct using batch_col
  expr.input.mb <- ComBat(dat = log2(expr.input.mb + 1), batch = clin.mb$RNA_library)
  expr.input.mb <- 2^expr.input.mb
}

# classify mb samples
print("Classify medulloblastoma subtypes...")
mb.classify <- classify.mat(expr.mb = expr.input.mb, clin.mb = clin.mb)

# save output to RData object
print("Writing output to file..")
saveRDS(mb.classify, file = output.file)
print("Done!")
