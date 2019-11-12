# Author: Komal S. Rathi
# Date: 10/07/2019
# Function: 
# 1. summarize RNA-seq to HUGO symbol x Sample matrix
# 2. tabulate corresponding gene annotations

# Example run: PolyA RNA-seq
# Rscript 01-summarize_matrices.R \
# -i ~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm.polya.rds \
# -g ~/Projects/OpenPBTA-analysis/data/raw/gencode.v27.primary_assembly.annotation.gtf \
# -m ~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds \
# -t ~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm-collapsed-table.polya.rds

# Example run: Stranded RNA-seq
# Rscript 01-summarize_matrices.R \
# -i ~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm.stranded.rds \
# -g ~/Projects/OpenPBTA-analysis/data/raw/gencode.v27.primary_assembly.annotation.gtf \
# -m ~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
# -t ~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm-collapsed-table.stranded.rds

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))

option_list <- list(
  make_option(c("-i", "--inputmat"), type = "character",
              help = "Input matrix of merged RSEM files (.RDS)"),
  make_option(c("-g", "--inputgtf"), type  = "character",
              help = "Input gtf file (.gtf)"),
  make_option(c("-m", "--outputmat"), type = "character",
              help = "Output matrix (.RDS)"),
  make_option(c("-t", "--dupstable"), type = "character",
              help = "Output gene annotation table (.RDS)")
)

# parse the parameters
opt <- parse_args(OptionParser(option_list = option_list))
input.mat <- opt$inputmat
input.gtf <- opt$inputgtf
output.mat <- opt$outputmat
dups.tab <- opt$dupstable

if(file.exists(output.mat) & file.exists(dups.tab)) {
  print("Matrix and Dups table exists!")
} else {
  print("Generating input matrix...!")
  expr <- readRDS(input.mat)
  
  # read gencode (v27 for the paper) to get gene id-gene symbo-gene type annotation
  # PAR_Y entries are just duplications - so remove beforehand
  gtf <- rtracklayer::import(input.gtf)
  gtf <- gtf %>% 
    as.data.frame() %>%
    select(gene_id, gene_name, gene_type) %>%
    mutate(gene_id = str_replace_all(gene_id, "_PAR_Y", "")) %>%
    unique()
  
  # split gene id and symbol
  expr <- expr %>% 
    mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
    separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
    unique() 
  
  # remove all genes with no expression
  expr <- expr[which(rowSums(expr[,3:ncol(expr)]) > 0),] 
  gtf$expressed <- ifelse(gtf$gene_id %in% expr$gene_id, "Yes", "No")
  
  # identify gene symbols mapped to multiple ensembl identifiers
  dups <- expr %>% 
    select(gene_id, gene_symbol) %>%
    group_by(gene_symbol) %>%
    mutate(gene_symbol.ct = n()) %>%
    filter(gene_symbol.ct > 1) %>%
    unique()
  gtf$ensembl_id <- ifelse(gtf$gene_id %in% dups$gene_id, "Multiple", "Unique")
  
  # collapse to matrix of HUGO symbols x Sample identifiers
  # take mean per row and use the max value for duplicated gene symbols
  expr.collapsed <- expr %>% 
    mutate(means = rowMeans(select(.,-gene_id, -gene_symbol))) %>% # take rowMeans
    arrange(desc(means)) %>% # arrange decreasing by means
    distinct(gene_symbol, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
    select(-means) %>%
    unique() %>%
    remove_rownames() 
  gtf$keep <- ifelse(gtf$gene_id %in% expr.collapsed$gene_id, "Yes", "No")
  
  # matrix of HUGO symbols x Sample identifiers
  expr.input <- expr.collapsed %>% 
    column_to_rownames("gene_symbol") %>%
    select(-c(gene_id)) 
  
  print(dim(expr.input))
  
  # final geneid-symbol-biotype annotation file
  print("Generating duplicates table...")
  saveRDS(object = gtf, file = dups.tab)
  
  # save matrix
  print("Saving collapsed matrix...")
  saveRDS(object = expr.input, file = output.mat)
  print("Done!!")
}
