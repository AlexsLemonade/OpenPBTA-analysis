# Author: Komal S. Rathi
# Date: 10/07/2019
# Function: 
# Script to summarize RNA-seq to HUGO symbol x Sample matrix
# Required input for immunedeconv

# Example run
# Rscript analyses/collapse-rnaseq/01-summarize_matrices.R -e data/pbta-gene-expression-rsem-fpkm.stranded.rds -o data/pbta-gene-expression-collapsed-rsem-fpkm.stranded.rds
# Rscript analyses/collapse-rnaseq/01-summarize_matrices.R -e data/pbta-gene-expression-rsem-fpkm.polya.rds -o data/pbta-gene-expression-collapsed-rsem-fpkm.polya.rds

# libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))

# input parameters
option_list <- list(
  make_option(c("-e", "--exprs"), type = "character",
              help = "Expression data (.RDS)"),
  make_option(c("-o","--outputmat"), type = "character",
              help = "Output matrix (.RDS)")
)

# parse the parameters
opt <- parse_args(OptionParser(option_list = option_list))
input.dat <- opt$expr
output.mat <- opt$outputmat
if(file.exists(output.mat)) {
  print("Matrix exists!")
} else {
  print("Generating input matrix...!")
  expr <- readRDS(input.dat)
  
  # reduce dataframe
  expr <- expr[which(rowSums(expr[,2:ncol(expr)]) > 0),] # remove all genes with no expression
  expr <- expr[grep('_PAR_', expr$gene_id, invert = T),] # discard PAR_* chromosomes from analysis
  
  # collapse to matrix of HUGO symbols x Sample identifiers
  # take mean per row and use the max value for duplicated gene symbols
  expr.collapsed <- expr %>% 
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    pivot_longer(-c(gene_id, gene_symbol), 
                 names_to = "sample_name", values_to = "fpkm") %>% 
    group_by(gene_id) %>%
    mutate(means = mean(fpkm)) %>% 
    group_by(gene_symbol) %>%
    filter(means == max(means)) %>% select(-gene_id) %>% unique()
  
  # matrix of HUGO symbols x Sample identifiers
  expr.input <- acast(expr.collapsed, gene_symbol~sample_name, value.var = 'fpkm')
  print(dim(expr.input))
  
  # save matrix
  saveRDS(object = expr.input, file = output.mat)
  print("Matrix generated. Done!!")
}
