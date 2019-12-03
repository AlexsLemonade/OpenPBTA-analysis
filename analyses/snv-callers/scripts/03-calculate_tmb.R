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
# --bed_wgs : File path that specifies the caller-specific BED regions file.
#             Assumes from top directory, 'OpenPBTA-analysis'.
# --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
#             is given from top directory of 'OpenPBTA-analysis'
# --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
#
# Command line example:
#
# Rscript analyses/snv-callers/scripts/03-calculate_tmb.R \
#   --consensus analyses/snv-callers/results/consensus/consensus_snv.maf.tsv \
#   --output analyses/snv-callers/results/consensus \
#   --metadata data/pbta-histologies.tsv \
#   --bed_wgs data/WGS.hg38.mutect2.unpadded.bed \
#   --bed_wxs data/WXS.hg38.100bp_padded.bed \
#   --no_region

################################ Initial Set Up ################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "wrangle_functions.R"))
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
    opt_str = "--bed_wgs", type = "character", default = "none",
    help = "File path that specifies the caller-specific
                BED regions file. Assumes from top directory, 'OpenPBTA-analysis'",
    metavar = "character"
  ),
  make_option(
    opt_str = "--bed_wxs", type = "character", default = "none",
    help = "File path that specifies the WXS BED regions file. Assumes
              from top directory, 'OpenPBTA-analysis'",
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
opt$db_file <- file.path(root_dir, opt$db_file)
opt$metadata <- file.path(root_dir, opt$metadata)
opt$bed_wgs <- file.path(root_dir, opt$bed_wgs)
opt$bed_wxs <- file.path(root_dir, opt$bed_wxs)

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$consensus, opt$metadata, opt$bed_wgs, opt$bed_wxs, opt$db_file 
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
tmb_coding_file <- file.path(opt$output, "consensus_snv_tmb_coding_only.tsv")
tmb_all_file <- file.path(opt$output, "consensus_snv_tmb_all.tsv")

########################### Set up this caller's data ##########################
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
# If the file exists or the overwrite option is not being used, run TMB calculations
if (file.exists(tmb_coding_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The Tumor Mutation Burden file already exists: \n",
    tmb_coding_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Print out warning if this file is going to be overwritten
  if (file.exists(tmb_coding_file)) {
    warning("Overwriting existing TMB file.")
  }
  # Print out progress message
  message(paste("Calculating TMB for", opt$consensus, "MAF data..."))
  
  # Set up BED region files for TMB calculations
  wgs_bed <- readr::read_tsv(opt$bed_wgs, col_names = FALSE)
  wxs_bed <- readr::read_tsv(opt$bed_wxs, col_names = FALSE)
  
  # Calculate size of genome surveyed
  wgs_genome_size <- sum(wgs_bed[, 3] - wgs_bed[, 2])
  wxs_exome_size <- sum(wxs_bed[, 3] - wxs_bed[, 2])
  
  # Print out these genome sizes
  cat(
    " WGS size in bp:", wgs_genome_size,
    "\n",
    "WXS size in bp:", wxs_exome_size,
    "\n"
  )
  
  # Only do this step if you have WXS samples
  if (any(metadata$experimental_strategy == "WXS")) {
    # Filter out mutations for WXS that are outside of these BED regions.
    maf_df <- wxs_bed_filter(maf_df, wxs_bed_file = opt$bed_wxs)
  }
  ############################ Coding TMB file #################################
  # Make a maf_df of only the coding mutations
  coding_maf_df <- maf_df %>%
    dplyr::filter(!(Variant_Classification %in% c("IGR", "Silent")))
  
  # Calculate coding only TMBs and write to file
  tmb_coding_df <- calculate_tmb(coding_maf_df,
                                 wgs_size = wgs_genome_size,
                                 wxs_size = wxs_exome_size
  )
  readr::write_tsv(tmb_coding_df, tmb_coding_file)
  
  # Print out completion message
  message(paste("TMB 'coding only' calculations saved to:", tmb_coding_file))
  
  ######################### All mutations TMB file #############################
  # Start up connection
  con <- DBI::dbConnect(RSQLite::SQLite(), opt$db_file)
  
  # Designate caller tables from SQL file
  strelka <- dplyr::tbl(con, "strelka")
  mutect <- dplyr::tbl(con, "mutect")
  
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
  strelka_mutect_maf_df <- strelka_mutect_maf_df  %>%
    dplyr::union_all(strelka_mutect_mnv, 
                     copy = TRUE) %>%
    dplyr::inner_join(metadata %>%
                        dplyr::select(
                          "Tumor_Sample_Barcode",
                          "experimental_strategy",
                          "short_histology"
                        ), 
                      copy = TRUE) %>% 
    as.data.frame()
  
  # Calculate TMBs and write to TMB file
  tmb_all_df <- calculate_tmb(strelka_mutect_maf_df,
                              wgs_size = wgs_genome_size,
                              wxs_size = wxs_exome_size
  )
  readr::write_tsv(tmb_all_df, tmb_all_file)
  
  # Print out completion message
  message(paste("TMB 'all' calculations saved to:", tmb_all_file))

}
