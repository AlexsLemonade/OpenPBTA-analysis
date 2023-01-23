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
`%>%` <- dplyr::`%>%`


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
  
  # tumor purity metadata file which has been filtered to only biospecimens that
  #  survive the cancer-group-level threshold
  tumor_purity_file <- file.path(
    rprojroot::find_root(rprojroot::has_dir(".git")), 
    "analyses", 
    "tumor-purity-exploration",
    "results",
    "thresholded_rna_stranded_same-extraction.tsv"
  )
  
  # Define vector of IDs to retain
  bs_ids <- readr::read_tsv(tumor_purity_file) %>%
    dplyr::pull(Kids_First_Biospecimen_ID) %>%
    unique()
  
  # Subset the `PTBA_GE_Standard` variable to those IDs
  PTBA_GE_Standard <- PTBA_GE_Standard[,bs_ids]
  
  # Quick check on filtering:
  if (!length(bs_ids) == ncol(PTBA_GE_Standard)) {
    stop("Error filtering IDs to use in EXTEND calculation with tumor purity thresholding on.")
  }
}


data = as.matrix(PTBA_GE_Standard)
RunEXTEND(data)#####EXTEND 
file.rename((pattern="TelomeraseScores.txt"), output_file)