# 00-subset-for-EPN.R
#
# Josh Shapiro for CCDL 2020
# 
# Purpose: Subsetting Expression data for EPN subtyping
#
# Option descriptions
# -h, --histology : path to the histology metadata file 
# -e, --expression : path to expression data file in RDS
# -o, --output_file : path for output file
#
# example invocation:
# Rscript scripts/bed_to_segfile.R \
#   -i results/cnv_consensus.tsv \
#   -o results/pbta-cnv-consensus.seg
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
    help = "hisology file tsv",
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
    help = "output file"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# read in files
histologies <- readr::read_tsv(opts$histology, col_types = readr::cols(.default = "c"))
expression <- readr::read_rds(opts$expression)


epr_samples <- histologies %>%
  filter(experimental_strategy == "RNA-Seq",
         disease_type_new == "Ependymoma") %>%
  pull(Kids_First_Biospecimen_ID)

# Subsetting expression columns with column names/BSIDs that are  in the  list  of ependymoma samples
epr_expression <- expression %>%
  select(epr_samples) %>%
  tibble::rownames_to_column("GENE")

# write the expression file out
readr::write_tsv(epr_expression, opts$outfile)
   

  
  