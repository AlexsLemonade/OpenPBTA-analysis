################################################################################
# This script conducts gene set enrichment analysis, specifically using the GSVA method [1] for scoring hallmark human pathway enrichment from RNA-Seq results. 
#
# The GSVA scores (i.e., enrichment scores) are calculated to produce a **Gaussian distribution of scores** "under the null hypothesis of no change in the pathway activity throughout the sample population."
# The authors claim a benefit to this approach:
#   + "Penalizes deviations that are large in *both* tails"
#   + "Provides a 'normalization' of the enrichment score by subtracting potential noise
#   + "Emphasizes genes in pathways that are concordantly activated in one direction only"
#   + "For pathways containing genes strongly acting in both directions, the deviations with cancel each other out and show little or no enrichment."
#
# Written by Stephanie J. Spielman for CCDL ALSF, 2020
#
#
#
# ####### USAGE, assumed to be run from top-level of project:
# Rscript --vanilla 'analyses/gene-set-enrichment-analysis/01-conduct-gsea-analysis.R --input <expression input file> --output <output file for writing scores>
#     --input : The name of the input expression data file to use for calculating scores. This file is assumed to live in the project's `data/` directory.
#     --output: The name of the TSV-formatted output file of GSVA scores, to be saved in the `results/` directory of this analysis.
# 
# Reference:
# 1. Sonja Hänzelmann, Robert Castelo, and Justin Guinney. 2013. “GSVA: Gene Set Variation Analysis for Microarray and RNA-Seq Data.” BMC Bioinformatics 14 (1): 7. https://doi.org/10.1186/1471-2105-14-7.
################################################################################



#### Set Up Libraries --------------------------------------------------------------------

## Load and/or install libraries ##
library(tidyr)
library(readr)
library(tibble)
library(optparse)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

if (!("msigdbr" %in% installed.packages())) {
  install.packages("msigdbr")
}
if (!("BiocManager" %in% installed.packages())) {
  install.packages("BiocManager")
}
library(BiocManager, quietly = TRUE)
if (!("GSVA" %in% installed.packages())) {
  BiocManager::install("GSVA")
}
library(msigdbr) ## Contains the hallmark data sets
library(GSVA)    ## Performs GSEA analysis


#### Set Up optparse --------------------------------------------------------------------

## Define arguments
option_list <- list(
  optparse::make_option(
    c("--input"),
    type = "character",
    default = NA,
    help = "The input file of expression data from which scores will be calculated."
  ),
    optparse::make_option(
    c("--output"),
    type = "character",
    default = NA, 
    help = "The output file for writing GSVA scores in TSV format."
  )
)

## Read in arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
if (is.na(opt$input)) stop("\n\nERROR: You must provide an input file with expression data with the flag --input, assumed to be in the `data/` directory of the OpenPBTA repository..") 
if (is.na(opt$output)) stop("\n\nERROR: You must provide an output file for saving GSVA scores with the flag --output, assumed to be placed in the `results/` directory of this analysis.") 


#### Set Up paths and file names --------------------------------------------------------------------

## Define directories
data_dir    <- file.path("..", "..", "data") 
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir)

## Ensure the input file exists in `data/` and specify input/output files 
expression_data_file <- file.path(data_dir, basename(opt$input))
if (!file.exists(expression_data_file)) stop("\n\nERROR: Provided input file does not exist.")
scores_output_file <- file.path(results_dir, basename(opt$output))



#### Load input files --------------------------------------------------------------------
expression_data <- as.data.frame( readr::read_rds(expression_data_file) )
human_hallmark  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") ## human hallmark genes from `migsdbr` package. The loaded data is a tibble.



#### Perform gene set enrichment analysis --------------------------------------------------------------------

# Prepare expression data: log2 transform re-cast as matrix
### Rownames are genes and column names are samples
expression_data_log2_matrix <- as.matrix( log2(expression_data + 1) )

# Prepare hallmark genes: Create a list of hallmarks, each of which is a list of genes
human_hallmark_twocols <- human_hallmark %>% dplyr::select(gs_name, human_gene_symbol)
human_hallmark_list    <- base::split(human_hallmark_twocols$human_gene_symbol, list(human_hallmark_twocols$gs_name))


#We then calculate the Gaussian-distributed scores
gsea_scores <- GSVA::gsva(expression_data_log2_matrix,
                          human_hallmark_list,
                          method = "gsva",
                          min.sz=1, max.sz=1500, ## Arguments from K. Rathi
                          mx.diff = TRUE)        ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)

### Clean scoring into tidy format
gsea_scores_df <- as.data.frame(gsea_scores) %>%
                    rownames_to_column(var = "hallmark_name")
  
#first/last_bs needed for use in gather (we are not on tidyr1.0)
first_bs <- head(colnames(gsea_scores), n=1)
last_bs  <- tail(colnames(gsea_scores), n=1)
  
gsea_scores_df_tidy <- gsea_scores_df %>%
                        tidyr::gather(Kids_First_Biospecimen_ID, gsea_score, !!first_bs : !!last_bs) %>%
                        dplyr::select(Kids_First_Biospecimen_ID, hallmark_name, gsea_score) 


#### Export GSEA scores to TSV --------------------------------------------------------------------
write_tsv(gsea_scores_df_tidy, scores_output_file)

