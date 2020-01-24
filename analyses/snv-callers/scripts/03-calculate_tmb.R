# Calculate TMB for a given MAF file
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# --consensus : File path to the MAF-like file.
# --db_file : Path to sqlite database file made from 01-setup_db.py
# --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
#              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --all_bed_wgs : File path that specifies the BED regions file to be used for the 
#                 denominator for all mutations TMB for WGS samples.
# --all_bed_wxs : File path that specifies the BED regions file to be used for the 
#                 denominator for all mutations TMB for WXS samples.
# --coding_bed_wgs : File path that specifies the BED regions file to be used for the 
#                 denominator for coding only TMB for WGS samples.
# --coding_bed_wxs : File path that specifies the BED regions file to be used for the 
#                 denominator for coding only TMB for WXS samples.
# --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
#
# Command line example:
#
# Rscript analyses/snv-callers/scripts/03-calculate_tmb.R \
# --consensus analyses/snv-callers/results/consensus/consensus_snv.maf.tsv \
# --db_file scratch/testing_snv_db.sqlite \
# --output analyses/snv-callers/results/consensus \
# --metadata data/pbta-histologies.tsv \
# --all_bed_wgs scratch/intersect_strelka_mutect_WGS.bed \
# --all_bed_wxs data/WXS.hg38.100bp_padded.bed \
# --coding_bed_wgs scratch/intersect_cds_lancet_strelka_mutect_WGS.bed \
# --coding_bed_wxs scratch/intersect_cds_WXS.bed \
# --overwrite

################################ Initial Set Up ################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "tmb_functions.R"))
source(file.path(root_dir, "analyses", "snv-callers", "util", "split_mnv.R"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-c", "--consensus"), type = "character",
    default = NULL, help = "File path to the MAF-like file",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-d", "--db_file"), type = "character",
    default = NULL, help = "Path to sqlite database file made from 01-setup_db.py",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL, help = "Path to folder where you would like the
              output TMB file from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--all_bed_wgs", type = "character", default = "none",
    help = "File path that specifies the BED regions file to be used for the 
    denominator for all mutations TMB for WGS samples.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--all_bed_wxs", type = "character", default = "none",
    help = "File path that specifies the BED regions file to be used for the 
    denominator for all mutations TMB for WXS samples.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--coding_bed_wgs", type = "character", default = "none",
    help = "File path that specifies the BED regions file to be used for the 
    denominator for coding only TMB for WXS samples. 'OpenPBTA-analysis'",
    metavar = "character"
  ),
  make_option(
    opt_str = "--coding_bed_wxs", type = "character", default = "none",
    help = "File path that specifies the BED regions file to be used for the 
    denominator for coding only TMB for WXS samples. 'OpenPBTA-analysis'",
    metavar = "character"
  ),
  make_option(
    opt_str = "--overwrite", action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any files of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Make everything relative to root path
opt$consensus <- file.path(root_dir, opt$consensus)
opt$metadata <- file.path(root_dir, opt$metadata)
opt$db_file <- file.path(root_dir, opt$db_file)
opt$all_bed_wgs <- file.path(root_dir, opt$all_bed_wgs)
opt$all_bed_wxs <- file.path(root_dir, opt$all_bed_wxs)
opt$coding_bed_wgs <- file.path(root_dir, opt$coding_bed_wgs)
opt$coding_bed_wxs <- file.path(root_dir, opt$coding_bed_wxs)

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$consensus, opt$metadata, opt$db_file, opt$all_bed_wgs, opt$all_bed_wxs, 
  opt$coding_bed_wgs, opt$coding_bed_wxs
)

# Get list of which files were found
files_found <- file.exists(needed_files)

# Report error if any of them aren't found
if (!all(files_found)) {
  stop(paste("\n Could not find needed file(s):",
    needed_files[which(!files_found)],
    "Check your options and set up.",
    sep = "\n"
  ))
}

############################### Set Up Output #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make output folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Declare output file paths
tmb_coding_file <- file.path(opt$output, "pbta-snv-consensus-mutation-tmb-coding.tsv")
tmb_all_file <- file.path(opt$output, "pbta-snv-consensus-mutation-tmb-all.tsv")

# Don't bother if both files exist already and overwrite is FALSE
if (all(file.exists(c(tmb_coding_file, tmb_all_file)), !opt$overwrite)) {
  stop(paste0(tmb_coding_file, tmb_all_file, "both exist and --overwrite is not
              being used. Use --overwrite if you would like to overwrite these files."))
}
########################### Set up this consensus data ##########################
# Print progress message
message(paste("Reading in", opt$consensus, "MAF data..."))

# Read in this MAF, skip the version number
maf_df <- data.table::fread(opt$consensus, data.table = FALSE)

# Print progress message
message("Setting up metadata...")

# Isolate metadata to only the samples that are in the datasets
metadata <- readr::read_tsv(opt$metadata) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% maf_df$Tumor_Sample_Barcode) %>%
  dplyr::distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE) %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID)

# Make sure that we have metadata for all these samples.
if (!all(unique(maf_df$Tumor_Sample_Barcode) %in% metadata$Tumor_Sample_Barcode)) {
  stop("There are samples in this MAF file that are not in the metadata.")
}

# Add the experimental strategy column on the data.frame for calculating purposes
maf_df <- maf_df %>%
  dplyr::inner_join(metadata %>%
    dplyr::select(
      Tumor_Sample_Barcode,
      experimental_strategy,
      short_histology
    ))

############################# Calculate TMB ####################################

############################# Coding TMB file ##################################
# If the file exists or the overwrite option is not being used, run TMB calculations
if (file.exists(tmb_coding_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The 'coding only' Tumor Mutation Burden file already exists: \n",
    tmb_coding_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Print out warning if this file is going to be overwritten
  if (file.exists(tmb_coding_file)) {
    warning("Overwriting existing 'coding only' TMB file.")
  }
  
  # Print out progress message
  message(paste("Calculating 'coding only' TMB..."))

  # Calculate coding only TMBs and write to file
  tmb_coding_df <- calculate_tmb(maf_df,
    bed_wgs = opt$coding_bed_wgs,
    bed_wxs = opt$coding_bed_wxs
  )
  readr::write_tsv(tmb_coding_df, tmb_coding_file)

  # Print out completion message
  message(paste("TMB 'coding only' calculations saved to:", tmb_coding_file))
}

########################### All mutations TMB file #############################
# If the file exists or the overwrite option is not being used, run TMB calculations
if (file.exists(tmb_all_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The 'all mutations' Tumor Mutation Burden file already exists: \n",
    tmb_all_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Print out warning if this file is going to be overwritten
  if (file.exists(tmb_coding_file)) {
    warning("Overwriting existing 'all mutations' TMB file.")
  }

  ######## Obtain Mutect Strelka mutations
  # Start up connection
  con <- DBI::dbConnect(RSQLite::SQLite(), opt$db_file)

  # Designate caller tables from SQL file
  strelka <- dplyr::tbl(con, "strelka")
  mutect <- dplyr::tbl(con, "mutect")

  # Specify the columns to join by
  join_cols <- c(
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Allele",
    "Tumor_Sample_Barcode"
  )

  # Create the consensus for non-MNVs
  strelka_mutect_maf_df <- strelka %>%
    dplyr::inner_join(mutect, by = join_cols)

  # Get Multi-nucleotide calls from mutect as SNVs
  split_mutect_df <- split_mnv(mutect)

  # join MNV calls with strelka
  strelka_mutect_mnv <- strelka %>%
    dplyr::inner_join(split_mutect_df,
      by = join_cols,
      copy = TRUE
    ) %>%
    as.data.frame()

  # Add in the MNVs
  strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
    dplyr::union_all(strelka_mutect_mnv,
      copy = TRUE
    ) %>%
    dplyr::inner_join(metadata %>%
      dplyr::select(
        "Tumor_Sample_Barcode",
        "experimental_strategy",
        "short_histology"
      ),
    copy = TRUE
    ) %>%
    as.data.frame()

  # .x is messing up the maf_to_granges function
  colnames(strelka_mutect_maf_df) <- gsub("\\.x$", "", colnames(strelka_mutect_maf_df))

  # Calculate TMBs and write to TMB file
  tmb_all_df <- calculate_tmb(strelka_mutect_maf_df,
    bed_wgs = opt$all_bed_wgs,
    bed_wxs = opt$all_bed_wxs
  )
  readr::write_tsv(tmb_all_df, tmb_all_file)

  # Print out completion message
  message(paste("TMB 'all' calculations saved to:", tmb_all_file))
}
