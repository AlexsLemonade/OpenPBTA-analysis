# K. S. Gaonkar 2021
# In this script, we will filter maf from a caller into 
# a resulting filtered file that has all calls that 
# - overlap with amino acid positions in a curated and 
# published cancer hotspot [database](https://www.cancerhotspots.org/files/hotspots_v2.xls)
# - overlap with non-coding region hotspots mutations ,


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))

option_list <- list(
  make_option(c("-m", "--maffile"),type="character",
              help="Filtered calls from [mutect2|strelka2|vardict|lancet]"),
  make_option(c("-c", "--caller"), type="character",
              help="Caller type [mutect2|strelka2|vardict|lancet]"),
  make_option(c("-f","--cancer_hotspot_folder"),
              help="MSKCC hotspot file", type="character"),
  make_option(c("-g","--genomic_site_hotspot_file"),
              help="genomic hotspot location file", type="character")
)

# Get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))
caller <- opt$caller
# Converting caller ID to call Uppercase for error handling
caller<-tolower(caller)

maf_coltypes <- readRDS(file.path("input","maf_coltypes.RDS"))
# maf_coltypes is a named vector with column types information for each column in a maf file
# More description about maf format is available 
# [here](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)

maf_df <- data.table::fread(opt$maffile,
                             stringsAsFactors = FALSE,
                             col.names = names(maf_coltypes),
                             colClasses = maf_coltypes)

mskcc_cancer_hotspot_file <- opt$cancer_hotspot_file
genomic_region_file <- opt$genomic_site_hotspot_file


# Directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
results_dir <- file.path(root_dir,
                         "analyses",
                         "hotspots-detection",
                         "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# MSKCC cancer hotspot file is divided into snv and indels
# in 2 tabs in the excel file provided in their [database](https://www.cancerhotspots.org/#/download) 
# We will combine the 2 sheets to create a `hotspot_database_amino_acid` dataframe 

# MSKCC cancer database
hotspot_database_2017_snv <- read_tsv(file.path(cancer_hotspot_folder,"hotspot_database_2017_snv.tsv")) 
hotspot_database_2017_indel <- read_tsv(file.path(cancer_hotspot_folder,"hotspot_database_2017_indel.tsv")) 
hotspot_database_amino_acid <- bind_rows( hotspot_database_2017_snv, hotspot_database_2017_indel) %>%
  dplyr::select("Amino_Acid_Position","Hugo_Symbol") %>%
  dplyr::mutate(hotspot_database="MSKCC") %>%
  unique() %>%
  as.data.frame()

# TERT promoter region
hotspot_database_genomic <- read_tsv(genomic_region_file)

# Gather Hugo_Symnol from both cancer hotspot and genomic regions to filter maf files
# This way we can directly filter large mafs into smaller manageable sizes to 
# check for overlaps
gene_table <- data.frame("Hugo_Symbol"=c(hotspot_database_amino_acid$Hugo_Symbol,
                                         hotspot_database_genomic$Hugo_Symbol))

print(paste("Subsetting" ,caller, "maf for hotspots"))

source ("utils/filterMaf.R")
# prepMaf will filter the input given maf file, given a
# - dataframe of hotspot_amino_acid_position_df to 
# check for Amino_Acid_Position overlaps per Hugo_Symbol
# - dataframe of hotspot_genomic_site_df to check for
# region overlap in the given Start_Position and End_Position
# - impact_values to filter IMPACT column in maf file
# - gene_table to filter Hugo_Symbol in maf file

maf_subset<-filterMaf(maf_df, gene_table = gene_table$Hugo_Symbol,
                        impact_values = c('HIGH', 'MODERATE', 'MODIFIER'),
                        hotspot_amino_acid_position_df = hotspot_database_amino_acid, 
                        hotspot_genomic_site_df = hotspot_database_genomic) %>%
  mutate(caller=caller)	

# save 
saveRDS(maf_subset, file.path(results_dir,paste0(caller,"_hotspots.RDS")))
