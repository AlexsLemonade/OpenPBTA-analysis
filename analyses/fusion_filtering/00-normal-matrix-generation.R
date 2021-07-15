# R. Jin 2021
# Filter expression matrix containing all specimens to only normal specimens of specific type

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))


option_list <- list(
  make_option(c("-e","--expressionMatrix"),type="character",
              help="expression matrix (TPM for samples that need to be zscore normalized .RDS)"),
  make_option(c("-c","--clinicalFile"),type="character",
              help="Histology file for all samples (.TSV) "),
  make_option(c("-s","--specimenType"),type="character",
              help="which specimen type to collect normal expression"),
  make_option(c("-o","--outputFile"),type="character",
              help="normalization TPM expression data to compare (.rds)")
)

opt <- parse_args(OptionParser(option_list=option_list))
expressionMatrix<-opt$expressionMatrix
clinicalFile<-opt$clinicalFile
specimenType<- opt$specimenType
outputFile<- opt$outputFile


#read in clinical file to find the list of normal specimens
histology_df <- read.delim(clinicalFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
normal_specimen <- histology_df %>% filter(gtex_group == specimenType) %>% 
  filter(experimental_strategy == "RNA-Seq") %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

#read in expression matrix with all specimens and only select 
expressionMatrix <- readRDS(expressionMatrix)
normal_expression_matrix <- expressionMatrix %>% select(rownames(normal_specimen))

saveRDS(normal_expression_matrix,outputFile)
