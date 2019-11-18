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
# --output_file : File path and file name of where you would like the MAF-like 
#                 output from this script to be stored.
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
# --vaf_filter 0.15 \
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
    opt_str = c("-o", "--output_file"), type = "character",
    default = NULL, help = "File path and file name of where you would like the 
                            MAF-like output from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--vaf_filter", type = "double", default = "0",
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

vaf_filter <- opt$vaf_filter # get out of opt list for sql

############################## Custom Functions ################################

#' Split multinucleotide variants into single nucleotide calls
#'
#' @param mnv_tbl a table containing MNVs (may be from an sql connection)
#'
#' @return a data frame of SNVs derived from the MNVs 
#'   (does not contain any of the original SNVs)
split_mnv <- function(mnv_tbl){
  mnv_df <- mnv_tbl %>%
    dplyr::filter(Variant_Type %in% c("DNP", "TNP", "ONP")) %>%
    as.data.frame() %>%
    # temp_id is not currently used, but will be handy for reconstituting MNVs
    dplyr::mutate(temp_id = dplyr::row_number()) %>%
    # lancet adds a base to the start of MNV Allele Calls
    dplyr::mutate(
      Allele = dplyr::case_when(
        stringr::str_length(Allele) == stringr::str_length(Reference_Allele) ~
          Allele,
        stringr::str_length(Allele) == stringr::str_length(Reference_Allele) + 1 ~ 
          stringr::str_sub(Allele, 2)
      )
    )
    # separate into individual rows
  if ("Match_Norm_Seq_Allele1" %in% colnames(mnv_df)){
    mnv_df <- mnv_df %>% 
      tidyr::separate_rows(Reference_Allele, 
                           Tumor_Seq_Allele1,
                           Tumor_Seq_Allele2,
                           Match_Norm_Seq_Allele1,
                           Match_Norm_Seq_Allele2,
                           Allele, 
                           sep = "(?<=[A-Za-z])")
  }else{
    mnv_df <- mnv_df %>% 
      tidyr::separate_rows(Reference_Allele, 
                           Tumor_Seq_Allele1,
                           Tumor_Seq_Allele2,
                           Allele, 
                           sep = "(?<=[A-Za-z])")
  }
  mnv_df <- mnv_df %>%
    # character separation leaves an extra blank
    dplyr::filter(Reference_Allele != "", 
                  Allele != "") %>%
    dplyr::group_by(temp_id) %>%
    dplyr::mutate(mnv_pos = dplyr::row_number(), 
                  Start_Position = Start_Position + mnv_pos - 1,
                  End_Position = Start_Position)
  return(mnv_df)
}


############################## Connect to database #############################
# Normalize this file path
opt$db_file <- file.path(root_dir, opt$db_file)

# Check that the database specified exists
if (!file.exists(opt$db_file)) {
  stop(paste("Error:", opt$db_file, "does not exist"))
}

# Start up connection
con <- DBI::dbConnect(RSQLite::SQLite(), opt$db_file)

############################### Set Up Output #####################################
# Normalize file path
opt$output_file <- file.path(root_dir, opt$output_file)

# Make sure the folder is made
output_dir <- stringr::word(opt$output_file, sep = "/", start = 1, end = -2)

# Make output folder
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

########################### Check for output file ###############################
# If the file exists or the overwrite option is not being used, do not write the
# merged VAF file.
if (file.exists(opt$output_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  stop(cat(
    "The merged consensus file already exists in the output directory: \n",
    opt$output_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
}

# Designate caller tables from SQL file
strelka <- dplyr::tbl(con, "strelka")
lancet <- dplyr::tbl(con, "lancet")
mutect <- dplyr::tbl(con, "mutect")

# We won't use VarDicts calls for the consensus
# vardict <- dplyr::tbl(con, "vardict")

# Use this so we can make columns back to traditional MAF order after join
full_col_order <- colnames(strelka)

# Specify the columns to join by
join_cols <- c(
  "Chromosome",
  "Start_Position",
  "Reference_Allele",
  "Allele",
  "Tumor_Sample_Barcode"
)

# Create the consensus for non-MNVs
consensus_df <- mutect %>%
  dplyr::inner_join(lancet, by = join_cols) %>%
  dplyr::select(join_cols) %>%
  dplyr::inner_join(strelka, by = join_cols) %>% 
  dplyr::select(full_col_order, -index)

# Get Multi-nucleotide calls from mutect and lancet as SNVs
split_mutect_df <- split_mnv(mutect)
split_lancet_df <- split_mnv(lancet)  

# join MNV calls with strelka
consensus_mnv <- strelka %>%
  dplyr::inner_join(split_mutect_df, 
                    by = join_cols,
                    copy = TRUE,
                    suffix = c("", "_mutect")) %>%
  dplyr::inner_join(split_lancet_df, 
                    by = join_cols,
                    copy = TRUE,
                    suffix = c("", "_lancet")) %>%
  dplyr::select(full_col_order, -index)

# MNV calls are not currently reconstituted from SNV representation, though this
# could be done at this point.


# Write consensus to output file
consensus_df %>%
  dplyr::union_all(consensus_mnv) %>%
  dplyr::arrange("Chromosome", "Start_Position") %>%
  # Filter with --vaf_filter option 
  dplyr::filter(VAF >= vaf_filter) %>% 
  as.data.frame() %>%
  # Write to a TSV file, change NAs back to "."
  readr::write_tsv(opt$output_file, na = ".")
