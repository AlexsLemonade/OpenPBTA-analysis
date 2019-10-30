# Merge the caller VAF and TMB files
#
# 2019
#
# C. Savonen for ALSF - CCDL
#
# Purpose: Merge callers' TMB and VAF files into total files with a column `caller`
# to designate their origin.

# Files Output:
# "all_callers_vaf.<file_format>" - contains all the VAF file information for all callers.
# "all_callers_tmb.<file_format>" - contains all the TMB file information for all callers.
# "mutation_id_list.<file_format>" - a full list of the mutations that can be
#                                    used for an UpSetR graph
# "callers_per_mutation.<file_format>" - contains a breakdown for each mutation of what callers
#                                        called it. Will be used to identify the consensus mutations.

# Option descriptions
# --vaf : Parent folder containing the vaf and tmb files for each folder.
#                                             <caller_name>_vaf.<file_format>
#                                             <caller_name>_tmb.<file_format>
# --file_format: What type of file format were the vaf and tmb files saved as? Options are
#               "rds" or "tsv". Default is "rds".
# --output : Where you would like the output from this script to be stored.
# --overwrite : If TRUE, will overwrite any reports of the same name. Default is
#              FALSE
#
#
# Command line example:
#
# Rscript 03-merge_callers.R \
# -v results \
# -o strelka2 \
# -s wxs \
# -w
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
    opt_str = c("--overwrite"), action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any reports of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

########################### Check options specified ############################
# Bring along the file suffix. Make to lower.
file_suffix <- tolower(opt$file_format)

# Check that the file format is supported
if (!(file_suffix %in% c("rds", "tsv"))) {
  warning("Option used for file format (-f) is not supported. Only 'tsv' or 'rds'
          files are supported. Defaulting to rds.")
  opt$file_format <- "rds"
  file_suffix <- "rds"
}

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
  value = TRUE
)

# Print this out to check
message("Will merge all VAF and TMB files in these folders: \n", paste0(caller_dirs, "\n"))

# Get a list of vaf files
vaf_files <- sapply(caller_dirs,
  list.files,
  pattern = paste0("_vaf.", file_suffix),
  recursive = TRUE, full.names = TRUE
)

# Print this out to check
message("Merging these VAF files: \n", paste0(vaf_files, "\n"))

# Get a list of tmb files
tmb_files <- sapply(caller_dirs,
  list.files,
  pattern = paste0("_tmb.", file_suffix),
  recursive = TRUE, full.names = TRUE
)

# Print this out to check
message("Merging these TMB files: \n", paste0(tmb_files, "\n"))

################################### Set Up #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make caller specific plots folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

########################### Make Master VAF file ###############################
# Get the caller names
caller_names <- stringr::word(vaf_files, sep = "/", -2)

# Read in vaf files for all callers
if (opt$file_format == "tsv") {
  vaf_list <- lapply(vaf_files, readr::read_tsv)
} else {
  vaf_list <- lapply(vaf_files, readr::read_rds)
}

# Read in the other files to match the first
vaf_list <- lapply(vaf_list, function(df) {
  # Make it so it is more easily combined with the other files
  df %>%
    # Attempt to make numeric columns where that doesn' kick back an "NA"
    dplyr::mutate_at(dplyr::vars(which(!is.na(as.numeric(t(df[1, ]))))), as.numeric) %>%
    # Aliquot id sometimes contains letters and sometimes numbers across the callers
    dplyr::mutate(
      aliquot_id = as.character(aliquot_id),
      variant_qual = as.character(variant_qual)
    ) %>%
    # Turn these columns into characters because otherwise they cause trouble.
    dplyr::mutate_at(dplyr::vars(dplyr::contains("AF", ignore.case = FALSE)), as.character) %>%
    # Get rid of the few if any duplicate entries.
    dplyr::distinct(mutation_id, .keep_all = TRUE)
})

# Carry over the callers' names
names(vaf_list) <- caller_names

# Print progress message
message("Saving master VAF file to: \n", file.path(opt$output, "all_callers_vaf.rds"))

# Combine and save VAF file
vaf_df <- suppressWarnings(dplyr::bind_rows(vaf_list, .id = "caller")) %>%
  dplyr::mutate(caller = factor(caller)) %>%
  # Write to RDS file
  readr::write_rds(file.path(opt$output, "all_callers_vaf.rds"))

########################### Make Master TMB file ###############################
# Read in TMB files for all callers
if (opt$file_format == "tsv") {
  tmb_list <- lapply(tmb_files, readr::read_tsv)
} else {
  tmb_list <- lapply(tmb_files, readr::read_rds)
}

# Carry over the callers' names
names(tmb_list) <- caller_names

# Print progress message
message("Saving master TMB file to: \n", file.path(opt$output, "all_callers_tmb.rds"))

# Combine and save TMB file
tmb_df <- dplyr::bind_rows(tmb_list, .id = "caller") %>%
  dplyr::mutate(caller = factor(caller)) %>%
  readr::write_rds(file.path(opt$output, "all_callers_tmb.rds"))

############################# Make mutation id list ############################
mutation_id_list <- lapply(vaf_list, function(caller) caller$mutation_id)

# Print progress message
message("Saving: \n", file.path(opt$output, "mutation_id_list.rds"))

readr::write_rds(mutation_id_list, file.path(opt$output, "mutation_id_list.rds"))

############################# Callers per mutation df ##########################
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

# Print progress message
message("Saving: \n", file.path(opt$output, "callers_per_mutation.rds"))

# Join the median VAF and the callers that call that mutation into one data.frame
callers_per_mutation <- callers_per_mutation %>%
  dplyr::inner_join(vaf_med, by = "mutation_id") %>%
  # Make column names more sensible
  dplyr::rename(caller_combo = "..x", median_vaf = "..y") %>%
  readr::write_rds(file.path(opt$output, "callers_per_mutation.rds"))
