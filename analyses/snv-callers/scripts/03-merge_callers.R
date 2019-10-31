# Merge the caller VAF and TMB files
#
# 2019
#
# C. Savonen for ALSF - CCDL
#
# Purpose: Merge callers' TMB and VAF files into total files with a column `caller`
# to designate their origin.

# Files Output:
# "all_callers_vaf.rds" - contains all the VAF file information for all callers.
# "all_callers_tmb.rds" - contains all the TMB file information for all callers.
# "mutation_id_list.rds" - a full list of the mutations that can be
#                                    used for an UpSetR graph
# "callers_per_mutation.rds" - contains a breakdown for each mutation of what callers
#                                        called it. Will be used to identify the consensus mutations.

# Option descriptions
# --vaf : Parent folder containing the vaf and tmb files for each folder.
#                                             <caller_name>_vaf.<file_format>
#                                             <caller_name>_tmb.<file_format>
# --output : Where you would like the output from this script to be stored.
# --overwrite : If TRUE, will overwrite any reports of the same name. Default is
#              FALSE
#
#
# Command line example:
#
# Rscript 03-merge_callers.R \
# -v results \
# -o results/consensus \
# --overwrite
#
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Read in RDS/TSV read in function
source(file.path(root_dir, "analyses", "snv-callers", "util", "read_function.R"))

# Load library:
library(optparse)

#--------------------------------Set up options--------------------------------#
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-v", "--vaf"), type = "character",
    default = NULL, help = "Path to folder with the output files
              from 01-calculate_vaf_tmb. Should include the VAF and TMB files",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL, help = "Path to folder where you would like the
              output from this script to be stored.",
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
# Normalize this file path
opt$vaf <- file.path(root_dir, opt$vaf)

# Check that the input directory exists
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
  pattern = "_vaf.",
  recursive = TRUE, full.names = TRUE
)

# Print this out to check
message("Merging these VAF files: \n", paste0(vaf_files, "\n"))

# Get a list of tmb files
tmb_files <- sapply(caller_dirs,
  list.files,
  pattern = "_tmb.",
  recursive = TRUE, full.names = TRUE
)

# Print this out to check
message("Merging these TMB files: \n", paste0(tmb_files, "\n"))

################################### Set Up #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make output folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Declare output file paths
all_vaf_file <- file.path(opt$output, "all_callers_vaf.rds")
all_tmb_file <- file.path(opt$output, "all_callers_tmb.rds")
mut_id_file <- file.path(opt$output, "mutation_id_list.rds")
call_per_mut_file <- file.path(opt$output, "callers_per_mutation.rds")

##################### Check for files if overwrite is FALSE ####################
# If overwrite is set to FALSE, check if these exist before continuing
if (!opt$overwrite) {
  # Make a list of the output files
  output_files <- c(all_vaf_file, all_tmb_file, mut_id_file, call_per_mut_file)

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

########################### Make Master VAF file ###############################
# If the file exists or the overwrite option is not being used, do not write the
# merged VAF file.
if (file.exists(all_vaf_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The merged VAF file already exists: \n",
    all_vaf_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Get the caller names
  caller_names <- stringr::word(vaf_files, sep = "/", -2)

  # Read in vaf files for all callers
  vaf_list <- lapply(vaf_files, read_tsv_or_rds)

  # Get the column names
  vaf_list_cols <- lapply(vaf_list, colnames)

  # Get the intersection of these column names
  vaf_list_cols <- base::Reduce(intersect, vaf_list_cols)

  # Read in the other files to match the first
  vaf_list <- lapply(vaf_list, function(df) {
    # Make it so it is more easily combined with the other files
    df <- df %>%
      # Only keep the columns that all the callers have
      dplyr::select(!!vaf_list_cols)
    df %>%
      # Attempt to make numeric columns where that doesn' kick back an "NA"
      dplyr::mutate_at(dplyr::vars(which(!is.na(as.numeric(t(df[1, ]))))), as.numeric) %>%
      # For some callers, aliquot id has characters, so we need all of them to be converted to chr
      dplyr::mutate(aliquot_id = as.character(aliquot_id)) %>%
      # Turn these columns into characters because otherwise they cause trouble.
      dplyr::mutate_at(dplyr::vars(dplyr::contains("AF", ignore.case = FALSE)), as.character) %>%
      # Get rid of the few if any duplicate entries.
      dplyr::distinct(mutation_id, .keep_all = TRUE)
  })

  # Carry over the callers' names
  names(vaf_list) <- caller_names

  # Print progress message
  message("Saving master VAF file to: \n", all_vaf_file)

  # Combine and save VAF file
  # Here `suppressWarnings` is being used because some caller VAF files do not have
  # certain annotation columns and their `NA`s or empty strings need to be coerced
  # so that they can be combined with the other callers.
  vaf_df <- suppressWarnings(dplyr::bind_rows(vaf_list, .id = "caller")) %>%
    dplyr::mutate(caller = factor(caller)) %>%
    # Write to RDS file
    readr::write_rds(all_vaf_file)
}
########################### Make Master TMB file ###############################
# If the file exists or the overwrite option is not being used, do not write the
# merged TMB file.
if (file.exists(all_tmb_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The merged TMB file already exists: \n",
    all_tmb_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  tmb_list <- lapply(tmb_files, read_tsv_or_rds)

  # Carry over the callers' names
  names(tmb_list) <- caller_names

  # Print progress message
  message("Saving master TMB file to: \n", all_tmb_file)

  # Combine and save TMB file
  tmb_df <- dplyr::bind_rows(tmb_list, .id = "caller") %>%
    dplyr::mutate(caller = factor(caller)) %>%
    readr::write_rds(all_tmb_file)
}
############################# Make mutation id list ############################
# If the file exists or the overwrite option is not being used, do not write mutation id file.
if (file.exists(mut_id_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The mutation id list file already exists: \n",
    mut_id_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  mutation_id_list <- lapply(vaf_list, function(caller) caller$mutation_id)

  # Print progress message
  message("Saving: \n", mut_id_file)

  readr::write_rds(mutation_id_list, mut_id_file)
}
############################# Callers per mutation df ##########################
# If the file exists or the overwrite option is not being used, do not write the
# callers per mutation file.
if (file.exists(call_per_mut_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The mutation id list file already exists: \n",
    call_per_mut_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
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
  message("Saving: \n", call_per_mut_file)

  # Join the median VAF and the callers that call that mutation into one data.frame
  callers_per_mutation <- callers_per_mutation %>%
    dplyr::inner_join(vaf_med, by = "mutation_id") %>%
    # Make column names more sensible
    dplyr::rename(caller_combo = "..x", median_vaf = "..y") %>%
    readr::write_rds(call_per_mut_file)
}
