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
# Rscript --vanilla 'analyses/gene-set-enrichment-analysis/01-conduct-gsea-analysis.R
#    Two optional commandline arguments:
#     --smallset=<TRUE/FALSE> (Default FALSE). Set to `TRUE` for running in CI
#     --smallset_size=<a numeric> (Default 15). Set the number of samples to run for smallset analysis, again for CI 
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
library(BiocManager)
if (!("GSVA" %in% installed.packages())) {
  BiocManager::install("GSVA")
}
library(msigdbr) ## Contains the hallmark data sets
library(GSVA)    ## Performs GSEA analysis


#### Set Up optparse --------------------------------------------------------------------

## Define arguments
option_list <- list(
  optparse::make_option(
    c("--smallset"),
    type = "logical",
    default = FALSE,
    help = "Set as TRUE to limit the number of samples evaluated, largely for running in CI."
  ),
    optparse::make_option(
    c("--smallset_size"),
    type = "double",
    default = 15,
    help = "Size of small sample, only used if argument `smallset` is specified as TRUE."
  )
)

## Read in arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


## Sanity check the options. 
MIN_SAMPLE_SIZE <- 15 # At least 15 should be run for small samples

if (is.null(opt$smallset) | is.na(opt$smallset)){
    opt$smallset <- FALSE
}
if (opt$smallset == 0) opt$smallset <- FALSE
if (opt$smallset == 1) opt$smallset <- TRUE

if (!(typeof(opt$smallset) == "logical"))  {
  stop("ERROR: The optional argument --smallset must be a logical argument, TRUE (1) or FALSE (0).")
}  

if (opt$smallset == TRUE){
  if (!(typeof(opt$smallset_size) == "double"))  {
    stop("ERROR: The optional argument --smallset_size must be a number")
  }
  # Ensure a reasonable amount at least
  if (opt$smallset_size < MIN_SAMPLE_SIZE) {
    opt$smallset_size <- MIN_SAMPLE_SIZE
  } 
  ## Ceil the smallset value in case of decimal
  opt$smallset_size <- ceiling(opt$smallset_size)
}

#### Set Up paths and file names --------------------------------------------------------------------

## Define directories
data_dir    <- file.path("..", "..", "data") 
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir)

## Define input files
## Define updated expression matrix. Columns are biospecimen ids and values are RSEM FPKM for stranded samples collapsed to gene symbol (gene-level)
expression_data_file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds")

## Define output file
gsea_scores_output_file <- file.path(results_dir, "gsva_scores.tsv")




#### Load input files --------------------------------------------------------------------
expression_data <- as.data.frame( readr::read_rds(expression_data_file) )
human_hallmark  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") ## human hallmark genes from `migsdbr` package. The loaded data is a tibble.



#### Perform gene set enrichment analysis --------------------------------------------------------------------

# Subset expression data if specified
if (opt$smallset) {
  expression_data <- expression_data[1:opt$smallset_size]
}

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
write_tsv(gsea_scores_df_tidy, gsea_scores_output_file)

