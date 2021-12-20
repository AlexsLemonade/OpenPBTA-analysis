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
              help = "EXTEND output (.txt)")
)



# parse the command line arguments - this becomes a list
# where each element is named by the long flag option
opt <- parse_args(OptionParser(option_list = option_list))
print(opt$input)
input_file <- opt$input
output_file <- opt$output



PTBA_GE_Standard = readRDS(input_file)
data = as.matrix(PTBA_GE_Standard)
RunEXTEND(data)#####EXTEND 
file.rename((pattern="TelomeraseScores.txt"), output_file)