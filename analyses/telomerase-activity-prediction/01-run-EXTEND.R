suppressPackageStartupMessages({
  library(pracma)
  library(dplyr)
  library(ggpubr)
  library(circlize)
  library(stringr)
  library(EXTEND)
  library(optparse)
})
stringsAsFactors=FALSE


###################################################   Running Analysis on RSEM-FPKM data formats  #####################################################################

# set up the command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Collapsed RSEM file (.rds)"),
  make_option(c("-o", "--output"), type = "character",
              help = "EXTEND output (.txt)"),
  make_option("--apply_tumor_purity_threshold",
              action = "store_true",
              default = FALSE,
              help = "This flag turns on filtering by tumor purity.")
)



# parse the command line arguments - this becomes a list
# where each element is named by the long flag option
opt <- parse_args(OptionParser(option_list = option_list))
print(opt$input)
input_file <- opt$input
output_file <- opt$output


PTBA_GE_Standard = readRDS(input_file)

# Filter data if tumor purity threshold is turned on
if (opt$apply_tumor_purity_threshold) {

  
  # Define path to tumor purity module
  tumor_purity_dir <- file.path(
    rprojroot::find_root(rprojroot::has_dir(".git")),
    "analyses",
    "tumor-purity-exploration"
  )
  
  # Define path to metadata file which has been filtered to only biospecimens that
  #  survive the cancer-group-level threshold
  tumor_purity_file <- file.path(
    tumor_purity_dir,
    "results",
    "thresholded_rna_stranded_same-extraction.tsv"
  )
  
  # Load the function to filter IDs
  source(
    file.path(tumor_purity_dir, "util", "function_filter-by-threshold.R")
  )
              
  
  # Filter the expression data
  PTBA_GE_Standard <- filter_expression_by_tumor_purity(PTBA_GE_Standard, 
                                                        tumor_purity_file)
    
}


data = as.matrix(PTBA_GE_Standard)
RunEXTEND(data)#####EXTEND
file.rename((pattern="TelomeraseScores.txt"), output_file)
