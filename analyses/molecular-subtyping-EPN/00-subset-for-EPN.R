# 00-subset-for-EPN.R
#
# Josh Shapiro for CCDL 2020
# 
# Purpose: Subsetting Expression data for EPN subtyping
#
# Option descriptions
# -h, --histology : path to the histology metadata file 
# -e, --expression : path to expression data file in RDS
# -o, --output_file : path for output tsv file, optionally gzipped with .gz ending
#
# example invocation:
# Rscript 00-subset-for-EPN.R \
#   -h ../../data/pbta-histologies.tsv \
#   -e ../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
#   -o epn-subset/epn-pbta-gene-expression-rsem-fpkm-collapsed.stranded.tsv.gz
# 

# Libraries
library(dplyr)
library(optparse)

# Parse command line options

option_list <- list(
  make_option(
    c("-i", "--histology"),
    type = "character",
    default = NULL,
    help = "histology file tsv",
  ),
  make_option(
    c("-e", "--expression"),
    type = "character",
    default = NULL,
    help = "expression file in RDS format",
  ),
  make_option(
    c("-o", "--outfile"),
    type = "character",
    default = NULL,
    help = "output tsv file; .gz for gzipped output."
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# read in files
histologies <- readr::read_tsv(opts$histology, col_types = readr::cols(.default = "c"))
expression <- readr::read_rds(opts$expression)


epn_samples <- histologies %>%
  filter(experimental_strategy == "RNA-Seq",
         disease_type_new == "Ependymoma") %>%
  pull(Kids_First_Biospecimen_ID)

# Subsetting expression columns with column names/BSIDs that are in the list of ependymoma samples
epn_expression <- expression %>%
  select(epn_samples) %>%
  tibble::rownames_to_column("GENE")

# write the expression file out
readr::write_tsv(epn_expression, opts$outfile)
   

  
  
