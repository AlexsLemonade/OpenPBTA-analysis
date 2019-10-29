# Consensus mutation file creation
# 2019
# C. Savonen for ALSF - CCDL
#
# Purpose: Merge callers' TMB and VAF files into total files with a column `caller`
# to designate their origin

# Option descriptions
# --label : Label to be used for folder and all output. eg. 'strelka2'. Optional.
#      Default is 'maf'
# --vaf : Parent folder containing the vaf and tmb files for each folder.
#                                             <caller_name>_vaf.<file_format>
#                                             <caller_name>_tmb.<file_format>
# --file_format: What type of file format were the vaf and tmb files saved as? Options are
#               "rds" or "tsv". Default is "rds".
# --output : Where you would like the output from this script to be stored.
# --overwrite : If TRUE, will overwrite any reports of the same name. Default is
#              FALSE
# --no_region : If used, regional analysis will not be done.

#
# Command line example:
#
# Rscript 02-run_eval.R \
# -l strelka2 \
# -v
# -o strelka2 \
# -s wxs \
# -w
#
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "merge_functions.R"))

# Load library:
library(optparse)

#--------------------------------Set up options--------------------------------#
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--label"), type = "character",
    default = "maf", help = "Label to be used for folder and all
                output. eg. 'strelka2'. Optional. Default is 'maf'",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-v", "--vaf"), type = "character",
    default = NULL, help = "Path to folder with the output files
              from 01-calculate_vaf_tmb. Should include the VAF, TMB, and
              region TSV files",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-f", "--file_format"), type = "character", default = "rds",
    help = "What type of file format were the vaf and tmb files saved as?
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
    opt_str = c("-w", "--overwrite"), action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any reports of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Bring along the file suffix. Make to lower.
file_suffix <- tolower(opt$file_format)

# Check that the file format is supported
if (!(file_suffix %in% c('rds', 'tsv'))) {
  warning("Option used for file format (-f) is not supported. Only 'tsv' or 'rds'
          files are supported. Defaulting to rds.")
  opt$file_format <- "rds"
  file_suffix <- "rds"
}
########################### Check options specified ############################
# Normalize this file path
opt$vaf <- file.path(root_dir, opt$vaf)

# Check the output directory exists
if (!dir.exists(opt$vaf)) {
  stop(paste("Error:", opt$vaf, "does not exist"))
}

# Exclude the non-caller directories
caller_dirs <- grep("vaf_cutoff|consensus",
                    dir(opt$vaf, full.names = TRUE),
                    invert = TRUE,
                    value = TRUE)

# Print this out to check
message("These are the callers whose files will be merged:", caller_dirs)

# The list of needed file suffixes
needed_files <- c(paste0(caller_dirs, "_vaf.", file_suffix),
                  paste0(caller_dirs,"_tmb.", file_suffix))

# Print the files we need
message(paste("Looking for these files:", needed_files))

# Get a list of vaf files
vaf_files <- sapply(caller_dirs,
                    list.files, pattern = paste0("_vaf.", file_suffix),
                     recursive = TRUE, full.names = TRUE)

# Get a list of tmb files
tmb_files <- sapply(caller_dirs,
                    list.files, pattern = paste0("_tmb.", file_suffix),
                     recursive = TRUE, full.names = TRUE)

# Report error if any of them aren't found
if(length(tmb_files) < length(caller_dirs)) {
  empty_dir <- caller_dirs[]
  stop(paste0(empty_dir, "does not have a ", "_tmb.", file_suffix, " file."))
}

# Specify the exact paths of these files
file_list <- file.path(opt$vaf, files_found)

################################### Set Up #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make caller specific plots folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

################################################################################
# Get the caller names
caller_names <- stringr::word(vaf_files, sep = "/", 2)

# Read in the other files to match the first
vaf_list <- lapply(vaf_files, function(vaf_file) {
  # Read in the file
  df <- readr::read_rds(vaf_file)

  # Make it so it is more easily combined with the other files
  df %>%
    # Attempt to make numeric columns where that doesn' kick back an "NA"
    dplyr::mutate_at(dplyr::vars(which(!is.na(as.numeric(t(df[1,]))))), as.numeric) %>%
    # Aliquot id sometimes contains letters and sometimes numbers across the callers
    dplyr::mutate(aliquot_id = as.character(aliquot_id),
                  variant_qual = as.character(variant_qual)) %>%
    # Turn these columns into characters because otherwise they cause trouble.
    dplyr::mutate_at(dplyr::vars(dplyr::contains("AF", ignore.case = FALSE)), as.character) %>%
    # Get rid of the few if any duplicate entries.
    dplyr::distinct(mutation_id, .keep_all = TRUE)
})

# Carry over the callers' names
names(vaf_list) <- caller_names

# Make master VAF
vaf_df <- dplyr::bind_rows(vaf_list, .id = "caller") %>%
  dplyr::mutate(caller = factor(caller)) %>%
  readr::write_rds(file.path(results_dir, "all_callers_vaf.rds"))

# Read in TMB files for all callers
tmb_list <- lapply(tmb_files, readr::read_rds)

# Carry over the callers' names
names(tmb_list) <- caller_names

tmb_df <- dplyr::bind_rows(tmb_list, .id = "caller") %>%
  dplyr::mutate(caller = factor(caller)) %>%
  readr::write_rds(file.path(results_dir, "all_callers_tmb.rds"))

# Make mutation id list
mutation_id_list <- lapply(vaf_list, function(caller) caller$mutation_id)
readr::write_rds(mutation_id_list, file.path(results_dir, "mutation_id_list.rds"))

# Make a string that says what callers call each mutation
callers_per_mutation <- tapply(vaf_df$caller,
  vaf_df$mutation_id,
  paste0,
  collapse = "-"
) %>%
  # Make into a data.frame
  as.data.frame() %>%
  tibble::rownames_to_column("mutation_id")

# Obtain the median VAF for each mutation
vaf_med <- tapply(
  vaf_df$vaf,
  vaf_df$mutation_id,
  median
) %>%
  # Make into a data.frame
  as.data.frame() %>%
  tibble::rownames_to_column("mutation_id")

# Join the median VAF and the callers that call that mutation into one data.frame
callers_per_mutation <- callers_per_mutation %>%
  dplyr::inner_join(vaf_med, by = "mutation_id")

# Make the column names more sensible
colnames(callers_per_mutation) <- c("mutation_id", "caller_combo", "median_vaf")
