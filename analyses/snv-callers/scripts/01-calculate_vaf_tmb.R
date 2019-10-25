# Run variant caller evaluation for a given MAF file.
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# -label : Label to be used for folder and all output. eg. 'strelka2'. Default is 'maf'.
# -output : File path that specifies the folder where the output should go.
#           New folder will be created if it doesn't exist. Assumes file path is
#           given from top directory of 'OpenPBTA-analysis'.
# --file_format: What type of file format would you like the output as? Options are
#               "rds" or "tsv". Default is "rds".
# --maf :  Relative file path to MAF file to be analyzed. Can be .gz compressed.
#          Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
#              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --annot_rds : Relative file path to annotation object RDS file to be analyzed.
#               Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --bed_wgs : File path that specifies the caller-specific BED regions file.
#             Assumes from top directory, 'OpenPBTA-analysis'.
# --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
#             is given from top directory of 'OpenPBTA-analysis'
# --vaf_filter: Optional Variant Allele Fraction filter. Specify a number; any
#               mutations with a VAF that are NA or below this number will be
#               removed from the vaf data.frame before it is saved to a TSV file.
# --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
# --no_region : If used, regional analysis will not be done. 
# --file_type
#
# Command line example:
#
# Rscript analyses/snv-callers/scripts/01-calculate_vaf_tmb.R \
#   --label strelka2 \
#   --output analyses/snv-callers/results \
#   --maf scratch/snv_dummy_data/strelka2 \
#   --metadata data/pbta-histologies.tsv \
#   --bed_wgs data/WGS.hg38.mutect2.unpadded.bed \
#   --bed_wxs data/WXS.hg38.100bp_padded.bed \
#   --annot_rds scratch/hg38_genomic_region_annotation.rds \
#   --vaf_filter .10 \
#   --no_region

################################ Initial Set Up ################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "wrangle_functions.R"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-l", "--label"), type = "character",
    default = "maf", help = "Label to be used for folder and all
                output. eg. 'strelka2'. Default is 'maf'",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character", default = "none",
    help = "File path that specifies the folder where the output should
              go. Assumes from top directory, 'OpenPBTA-analysis'. New folder
              will be created if it doesn't exist.",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-f", "--file_format"), type = "character", default = "rds",
    help = "What type of file format would you like the output as? Options are
            'rds' or 'tsv'. Default is 'rds'.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--maf", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-a", "--annot_rds"), type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to annotation object RDS file to be analyzed.",
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
    opt_str = "--vaf_filter", type = "character", default = "0",
    help = "Optional Variant Allele Fraction filter. Specify a number; any
            mutations with a VAF that are NA or below this number will be
            removed from the vaf data.frame before it is saved to a TSV file.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--overwrite", action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any files of
              the same name. Default is FALSE",
    metavar = "character"
  ),
  make_option(
    opt_str = "--no_region", action = "store_false",
    default = TRUE, help = "If used, regional analysis will not be run.",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Coerce to numeric
opt$vaf_filter <- as.numeric(opt$vaf_filter) 

# Make everything relative to root path
opt$maf <- file.path(root_dir, opt$maf)
opt$metadata <- file.path(root_dir, opt$metadata)
opt$bed_wgs <- file.path(root_dir, opt$bed_wgs)
opt$bed_wxs <- file.path(root_dir, opt$bed_wxs)
opt$cosmic <- file.path(root_dir, opt$cosmic)

# Bring along the file suffix. Make to lower. 
file_suffix <- tolower(opt$file_format)

# Check that the file format is supported
if (!(file_suffix %in% c('rds', 'tsv'))) {
  stop("Option used for file format (-f) is not supported. Only 'tsv' or 'rds'
       files are supported.")
}

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$maf, opt$metadata, opt$bed_wgs, opt$bed_wxs, opt$cosmic
)

# Only if regional analysis is being done do we need the annotation file
if (opt$no_region) {
  opt$annot_rds <- file.path(root_dir, opt$annot_rds)
  needed_files <- c(needed_files, opt$annot_rds)
}

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

################## Create output directories for this caller ##################
# Caller specific results directory path
caller_results_dir <- file.path(root_dir, opt$output)

# Make caller specific results folder
if (!dir.exists(caller_results_dir)) {
  dir.create(caller_results_dir, recursive = TRUE)
}

####################### File paths for files we will create ####################
vaf_file <- file.path(caller_results_dir, paste0(opt$label, "_vaf.", file_suffix))
region_annot_file <- file.path(caller_results_dir, paste0(opt$label, "_region.", file_suffix))
tmb_file <- file.path(caller_results_dir, paste0(opt$label, "_tmb.tsv", file_suffix))

# Declare metadata file name for this caller
metadata_file <- file.path(
  caller_results_dir,
  paste0(opt$label, "_metadata_filtered", file_suffix)
)

##################### Check for files if overwrite is FALSE ####################
# If overwrite is set to FALSE, check if these exist before continuing
if (!opt$overwrite) {
  # Make a list of the output files
  output_files <- c(vaf_file, region_annot_file, tmb_file)

  # Find out which of these exist
  existing_files <- file.exists(output_files)

  # If all files exist; stop
  if (all(existing_files)) {
    stop(cat(
      "Stopping; --overwrite is not being used and all output files already exist: \n",
      vaf_file, "\n",
      region_annot_file, "\n",
      tmb_file
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

########################### Set up this caller's data ##########################
# Print progress message
message(paste("Reading in", opt$maf, "MAF data..."))

# Read in this MAF, skip the version number
maf_df <- data.table::fread(opt$maf, skip = 1, data.table = FALSE)

# Print progress message
message(paste("Setting up", opt$label, "metadata..."))

# Isolate metadata to only the samples that are in the datasets
metadata <- readr::read_tsv(opt$metadata) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% maf_df$Tumor_Sample_Barcode) %>%
  dplyr::distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE) %>%
  dplyr::arrange() %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) 

# Write the metadata file
if (opt$file_format == "rds") {
  metadata %>%
    readr::write_rds(metadata_file)
} else {
  metadata %>%
    readr::write_tsv(metadata_file)
}

# Print out completion message
message(paste("Filtered metadata file saved to: \n", metadata_file))

# Make sure that we have metadata for all these samples.
if (!all(unique(maf_df$Tumor_Sample_Barcode) %in% metadata$Tumor_Sample_Barcode)) {
  stop("There are samples in this MAF file that are not in the metadata.")
}

################## Calculate VAF and set up other variables ####################
# If the file exists or the overwrite option is not being used, calculate VAF
if (file.exists(vaf_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The VAF file already exists: \n",
    vaf_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Print out warning if this file is going to be overwritten
  if (file.exists(vaf_file)) {
    warning("Overwriting existing VAF file.")
  }
  # Print out progress message
  message(paste("Calculating, sampling, and merging VAF for", opt$label, "MAF data..."))

  # Use the premade function to calculate VAF this will also merge the metadata
  vaf_df <- set_up_maf(maf_df, metadata, opt$vaf_filter)

  message(paste(nrow(vaf_df), "mutations left after filter and merge"))

  # Write the vaf file
  if (opt$file_format == "rds") {
    vaf_df %>% 
      readr::write_rds(vaf_file)
  } else {
    vaf_df %>% 
      readr::write_tsv(vaf_file)
  }
  # Print out completion message
  message(paste("VAF calculations saved to: \n", vaf_file))
}
######################### Annotate genomic regions #############################
if (opt$no_region) {
  # If the file exists or the overwrite option is not being used, run regional annotation analysis
  if (file.exists(region_annot_file) && !opt$overwrite) {
    # Stop if this file exists and overwrite is set to FALSE
    warning(cat(
      "The regional annotation file already exists: \n",
      region_annot_file, "\n",
      "Use --overwrite if you want to overwrite it."
    ))
  } else {
    # Print out warning if this file is going to be overwritten
    if (file.exists(region_annot_file)) {
      warning("Overwriting existing regional annotation file.")
    }
    # Print out progress message
    message(paste("Annotating genomic regions for", opt$label, "MAF data..."))

    # Annotation genomic regions
    maf_annot <- annotr_maf(vaf_df, annotation_file = opt$annot_rds) %>%
      readr::write_tsv(region_annot_file)

    # Print out completion message
    message(paste("Genomic region annotations saved to:", region_annot_file))
  }
}
############################# Calculate TMB ####################################
# If the file exists or the overwrite option is not being used, run TMB calculations
if (file.exists(tmb_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The Tumor Mutation Burden file already exists: \n",
    tmb_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Print out warning if this file is going to be overwritten
  if (file.exists(tmb_file)) {
    warning("Overwriting existing TMB file.")
  }
  # Print out progress message
  message(paste("Calculating TMB for", opt$label, "MAF data..."))

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
    vaf_df <- wxs_bed_filter(vaf_df, wxs_bed_file = opt$bed_wxs)
  }

  # Calculate TMBs and write to TMB file
  tmb_df <- calculate_tmb(vaf_df,
    wgs_size = wgs_genome_size,
    wxs_size = wxs_exome_size
  ) 
  
  # Write the tmb file
  if (opt$file_format == "rds") {
    tmb_df %>%
      readr::write_tsv(tmb_file)
  } else {
    tmb_df %>%
      readr::write_tsv(tmb_file)
  }

  # Print out completion message
  message(paste("TMB calculations saved to:", tmb_file))
}
