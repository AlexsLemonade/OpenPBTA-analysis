# Consensus mutation file creation
#
# 2019
#
# C. Savonen for ALSF - CCDL
#
# Purpose: Save consensus mutation calls to a MAF-like file.
#
# Option descriptions
# --merged_dir : File path to where the 03-merged_callers.R output can be found.
#                  Files required: "all_callers_vaf.rds"
#                                  "all_callers_tmb.rds"
#                                  "callers_per_mutation.rds"
# --vaf : What VAF should be used when combining callers? Options are
#         one of the caller names."
# --combo : What combination of callers need to detect a mutation for it to be
#           considered real and placed in the consensus file? List the callers
#           that need to be considered in alphabetical order with '-'
#           in between. eg. 'lancet-mutect2-strelka2'
# --output : Where you would like the output from this script to be stored.
# --bed_wgs : File path that specifies the caller-specific BED regions file.
#             Assumes from top directory, 'OpenPBTA-analysis'.
# --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
#             is given from top directory of 'OpenPBTA-analysis'
# --overwrite : If TRUE, will overwrite any reports of the same name. Default is
#              FALSE
#
# Command line example:
# Rscript 04-create_consensus_mut_files.R \
# --merged_dir results/consensus \
# --combo lancet-mutect2-strelka2 \
# --output analyses/snv-callers/results/consensus \
# --vaf strelka2 \
# --bed_wgs data/WGS.hg38.strelka2.unpadded.bed \
# --bed_wxs data/WXS.hg38.100bp_padded.bed \
# --overwrite
#
#################################### Set Up ####################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import this for the recalculation of TMB for the consensus mutation calls
source(file.path(root_dir, "analyses", "snv-callers", "util", "wrangle_functions.R"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-m", "--merged_dir"), type = "character",
    default = "", help = "File path to where the merged files were sent in
                          03-merge_callers.R",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL, help = "Path to folder where you would like the
              output from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = c("--vaf"), type = "character",
    default = FALSE, help = "What VAF should be used when combining callers?
    Options one of the caller names.",
    metavar = "character"
  ),
  make_option(
    opt_str = c("--combo"), type = "character",
    default = FALSE, help = "What combination of callers need to detect a
    mutation for it to be considered real and placed in the consensus file?
    List the callers that need to be considered in alphabetical order with '-'
    in between. eg. 'lancet-mutect2-strelka2'",
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
    default = FALSE, help = "If TRUE, will overwrite any reports of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

########################### Check options specified ############################
# Normalize these file paths
opt$merged_dir <- file.path(root_dir, opt$merged_dir)
opt$output <- file.path(root_dir, opt$output)
opt$bed_wgs <- file.path(root_dir, opt$bed_wgs)
opt$bed_wxs <- file.path(root_dir, opt$bed_wxs)

# Check the output directory exists
if (!dir.exists(opt$merged_dir)) {
  stop(paste("Error:", opt$merged_dir, "does not exist"))
}
# Check for BED files
if (!file.exists(opt$bed_wgs)) {
  stop(paste("Error:", opt$bed_wgs, "does not exist"))
}
# Check for BED files
if (!file.exists(opt$bed_wxs)) {
  stop(paste("Error:", opt$bed_wxs, "does not exist"))
}

# The list of needed file names but without their suffix
needed_files <- file.path(opt$merged_dir, c(
  "all_callers_vaf.rds",
  "all_callers_tmb.rds",
  "callers_per_mutation.rds"
))

# Report error if any of them aren't found
if (any(!file.exists(needed_files))) {
  stop(
    "Error: the directory specified with --output, doesn't have the ",
    "necessary file(s): \n", paste(needed_files[!file.exists(needed_files)], "\n")
  )
}
################################### Set Up #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make output folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Declare output file paths
consensus_mut_file <- file.path(opt$output, "consensus_mutation.maf.tsv")
consensus_mut_zip_file <- file.path(opt$output, "consensus_mutation.maf.tsv.zip")
consensus_tmb_file <- file.path(opt$output, "consensus_mutation.tmb.tsv")

##################### Check for files if overwrite is FALSE ####################
# If overwrite is set to FALSE, check if these exist before continuing
if (!opt$overwrite) {
  # Make a list of the output files
  output_files <- c(consensus_mut_file, consensus_mut_zip_file, consensus_tmb_file)

  # Find out which of these exist
  existing_files <- file.exists(output_files)

  # If all files exist; stop
  if (all(existing_files)) {
    stop(cat(
      "Stopping; --overwrite is not being used and all output files already exist
      in the designated --output directory."
    ))
  }
  # If some files exist, print a warning:
  if (any(existing_files)) {
    warning(cat(
      "Some output files already exist and will not be overwritten unless you use --overwrite: \n",
      paste0(output_files[which(existing_files)], "\n")
    ))
  }
}
############################### Read in the files ##############################
# Read in merged VAF_df
vaf_df <- readr::read_rds(input_files[[1]])

# If the caller chosen for the VAF reported isn't in the column, stop.
if (!(opt$vaf %in% vaf_df$caller)) {
  stop("Caller chosen with --vaf does not exist in the 'callers column in the
  master VAF file.")
}
# Read in merged TMB file
tmb_df <- readr::read_rds(input_files[[2]])

# Read in caller per mutation file
callers_per_mutation <- readr::read_rds(input_files[[3]])

# If the file exists or the overwrite option is not being used, do not write the
# consensus mutation file.
if (file.exists(consensus_mut_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "A consensus mutation file already exists: \n",
    consensus_mut_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Let's get the mutation ids for combination of
  consen_mutations <- callers_per_mutation %>%
    dplyr::filter(caller_combo == dplyr::sym(opt$combo)) %>%
    dplyr::select(-caller_combo) %>%
    dplyr::pull(mutation_id)

  # Stop if no mutations are found
  if (length(consen_mutations) == 0) {
    stop("No mutations had this combination of callers call it. Double check that
       the combination you specified with --combo is spelled exactly right and
       the callers are in alphabetical order. ")
  }

  # Isolate the mutations to only these mutations, use the Strelka2 stats.
  consen_mutation <- vaf_df %>%
    dplyr::filter(
      caller == opt$vaf,
      mutation_id %in% consen_mutations
    ) %>%
    readr::write_tsv(consensus_mut_file)
}

# If the file exists or the overwrite option is not being used, do not zip the
# consensus mutation file.
if (file.exists(consensus_mut_zip_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "A zipped consensus mutation file already exists: \n",
    consensus_mut_zip_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Zip this file up.
  zip(consensus_mut_zip_file, consensus_mut_file)
}
############################### Re-calculate TMB ###############################
# If the file exists or the overwrite option is not being used, do not create
# the consensus TMB file
if (file.exists(consensus_tmb_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "A consensus TMB file already exists: \n",
    consensus_tmb_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Set up BED region files for TMB calculations
  wgs_bed <- readr::read_tsv(file.path(opt$bed_wgs), col_names = FALSE)
  wxs_bed <- readr::read_tsv(file.path(opt$bed_wxs), col_names = FALSE)

  # Calculate size of genome surveyed
  # These files are BED files where the third column is the End position and 
  # the second column is the Start position. 
  # So End - Start gives the size of each range. Sum the gives the total size in bp.
  wgs_genome_size <- sum(wgs_bed[, 3] - wgs_bed[, 2])
  wxs_exome_size <- sum(wxs_bed[, 3] - wxs_bed[, 2])

  # Calculate TMBs and write to TMB file
  tmb_df <- calculate_tmb(vaf_df,
    wgs_size = wgs_genome_size,
    wxs_size = wxs_exome_size
  ) %>%
    readr::write_tsv(file.path(opt$output, "consensus_mutation_tmb.tsv"))

  # Give message
  message("Consensus mutations and TMB re-calculations saved in: \n", opt$output)
}
