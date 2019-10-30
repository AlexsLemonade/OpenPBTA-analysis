# Consensus mutation file creation
#
# 2019
#
# C. Savonen for ALSF - CCDL
#
# Purpose: Save consensus mutation calls to a MAF-like file.
#
# Option descriptions
# --merged_files : File path to where the 03-merged_callers.R output can be found.
#                  Files required: "all_callers_vaf.<file_format>"
#                                  "all_callers_tmb.<file_format>"
#                                  "mutation_id_list.<file_format>"
# --vaf : What VAF should be used when combining callers? Options are 'median' or
#         one of the caller names."
# --combo : What combination of callers need to detect a mutation for it to be
#           considered real and placed in the consensus file? List the callers
#           that need to be considered in alphabetical order with '-'
#           in between. eg. 'lancet-mutect2-strelka2'
# --file_format: What type of file format were the vaf and tmb files saved as? Options are
#               "rds" or "tsv". Default is "rds".
# --output : Where you would like the output from this script to be stored.
# --bed_wgs : File path that specifies the caller-specific BED regions file.
#             Assumes from top directory, 'OpenPBTA-analysis'.
# --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
#             is given from top directory of 'OpenPBTA-analysis'
# --overwrite : If TRUE, will overwrite any reports of the same name. Default is
#              FALSE

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
    opt_str = c("-m", "--merged_files"), type = "character",
    default = "", help = "File path to where the merged files were sent in 
                          03-merge_callers.R",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-f", "--file_format"), type = "character", default = "rds",
    help = "What type of file format were the files saved as?
            Options are 'rds' or 'tsv'. Default is 'rds'.",
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
    Options are 'median' or one of the caller names.",
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
    opt_str = c("--overwrite"), action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any reports of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

opt$merged_files <- "analyses/snv-callers/results/consensus"
opt$combo <- "lancet-mutect2-strelka2"
opt$vaf <- "strelka2"
opt$output <- "analyses/snv-callers/results/consensus"
opt$overwrite <- TRUE
opt$bed_wgs <- "data/WGS.hg38.strelka2.unpadded.bed"
opt$bed_wxs <- "data/WXS.hg38.100bp_padded.bed"

########################### Check options specified ############################
# Normalize these file paths
opt$merged_files <- file.path(root_dir, opt$merged_files)
opt$output <- file.path(root_dir, opt$output)
opt$bed_wgs <- file.path(root_dir, opt$bed_wgs)
opt$bed_wxs <- file.path(root_dir, opt$bed_wxs)

# Check the output directory exists
if (!dir.exists(opt$merged_files)) {
  stop(paste("Error:", opt$merged_files, "does not exist"))
}

# The list of needed file suffixes
needed_files <- c(
  "all_callers_vaf.",
  "all_callers_tmb.",
  "mutation_id_list.",
  "callers_per_mutation."
)

needed_files <- c(
  file.path(root_dir, opt$merged_files, paste0(needed_files, opt$file_format)),
  opt$bed_wgs,
  opt$bed_wxs
)

# Get list of which files were found
files_found <- sapply(needed_files, file.exists)

# Report error if any of them aren't found
if (any(is.na(files_found))) {
  stop(paste0(
    "Error: the directory specified with --output, doesn't have the",
    "necessary file(s):", names(files_found)[which(!files_found)]
  ))
}

############################### Read in the files ##############################
vaf_df <- readr::read_rds(file.path(
  opt$merged_files,
  paste0("all_callers_vaf.", opt$file_format)
))
tmb_df <- readr::read_rds(file.path(
  opt$merged_files,
  paste0("all_callers_tmb.", opt$file_format)
))
callers_per_mutation <- readr::read_rds(file.path(
  opt$merged_files,
  paste0("callers_per_mutation.", opt$file_format)
))

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
    caller == "strelka2",
    mutation_id %in% consen_mutations
  ) %>%
  readr::write_tsv(file.path(opt$output, "consensus_mutation.maf.tsv"))

# Zip this file up.
zip(
  file.path(opt$output, "consensus_mutation.maf.tsv.zip"),
  file.path(opt$output, "consensus_mutation.maf.tsv")
)

############################### Re-calculate TMB ###############################
# Set up BED region files for TMB calculations
wgs_bed <- readr::read_tsv(file.path(opt$bed_wgs), col_names = FALSE)
wxs_bed <- readr::read_tsv(file.path(opt$bed_wxs), col_names = FALSE)

# Calculate size of genome surveyed
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
