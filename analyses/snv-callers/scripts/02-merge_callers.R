# Make the consensus SNV file
#
# 2019
#
# C. Savonen for ALSF - CCDL
#
# Purpose: Merge callers' data files into consensus MAF-like file
#
# Option descriptions
# 
# --db_file : Path to sqlite database file made from 01-setup_db.py
# --output : Where you would like the output from this script to be stored.
# --vaf_filter: Optional Variant Allele Fraction filter. Specify a number; any
#               mutations with a VAF that are NA or below this number will be
#               removed from the vaf data.frame before it is saved to a TSV file.
# --overwrite : If TRUE, will overwrite any reports of the same name. Default is
#              FALSE
#
#
# Command line example:
#
# Rscript 01-merge_callers.R \
# --db_file results \
# --output results/consensus \
# --overwrite
#
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)

#--------------------------------Set up options--------------------------------#
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-d", "--db_file"), type = "character",
    default = NULL, help = "Path to sqlite database file made from 01-setup_db.py",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL, help = "Path to folder where you would like the
              output MAF-like file from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--vaf_filter", type = "character", default = "0",
    help = "Optional Variant Allele Fraction filter. Specify a number; any
            mutations with a VAF that are NA or below this number will be
            removed from the vaf data.frame before it is saved to a TSV file.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--overwrite", action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any reports of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

############################## Connect to database #############################
# Normalize this file path
opt$db_file <- file.path(root_dir, opt$db_file)

# Check that the database specified exists
if (!file.exists(opt$db_file)) {
  stop(paste("Error:", opt$db_file, "does not exist"))
}

# Start up connection 
con <- DBI::dbConnect(RSQLite::SQLite(), opt$db_file)

###############################Set Up Output #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make output folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Declare output file paths
consensus_maf_file <- file.path(opt$output, "consensus_snv.maf.tsv")

########################### Check for output file ###############################
# If the file exists or the overwrite option is not being used, do not write the
# merged VAF file.
if (file.exists(consensus_maf_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  stop(cat(
    "The merged consensus file already exists in the output directory: \n",
    opt$output, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
}

# Designate caller tables from SQL file
strelka <- dplyr::tbl(con, "strelka")
lancet <- dplyr::tbl(con, "lancet")
mutect <- dplyr::tbl(con, "mutect")

# We won't use VarDicts calls for the consensus
# vardict <- dplyr::tbl(con, "vardict")
  
# Specify the columns to join by
join_cols = c("Chromosome",
              "Start_Position",
              "Reference_Allele",
              "Allele",
              "Tumor_Sample_Barcode")
  
# Create the consensus
consensus_df <- strelka %>%
  dplyr::inner_join(lancet,
             by = join_cols,
             suffix = c("s", "l")) %>%
  dplyr::inner_join(mutect,
             by = join_cols,
             suffix = c("s", "m"))

# Write consensus to output file
consensus_df %>% 
  as.data.frame() %>% 
  readr::write_tsv(consensus_maf_file)
