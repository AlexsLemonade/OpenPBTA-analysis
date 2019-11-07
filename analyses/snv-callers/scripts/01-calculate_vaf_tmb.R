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
# --file_save: What type of file format would you like the output as? Options are
#               "rds" or "tsv". Default is "rds".
# --maf :  Relative file path to MAF file to be analyzed. Can be .gz compressed.
#          Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --sql_file : File path to where the SQL file was saved in 00-set_up.R.
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
#
# Command line example:
#
# Rscript analyses/snv-callers/scripts/01-calculate_vaf_tmb.R \
#   --label strelka2 \
#   --output analyses/snv-callers/results \
#   --maf scratch/snv_dummy_data/strelka2 \
#   --sql_file maf.sqlite \
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
    opt_str = c("-f", "--file_save"), type = "character", default = "none",
    help = "If you would like the file saved separately besides the sql file, 
            specify what you would like the output as. Options are
            'rds' or 'tsv'. Default is to not save any file besides the sqlite
            file.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--maf", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--sql_file", type = "character", default = "none",
    help = "File path to where the SQL file was saved in 00-set_up.R",
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

opt$label <- "strelka2"
opt$output <- "analyses/snv-callers/results/strelka2"
opt$file_save <- "none"
opt$maf <-  "data/testing/release-v9-20191105/pbta-snv-strelka2.vep.maf.gz"
opt$sql_file <- "analyses/snv-callers/ref_files/maf.sqlite"
opt$bed_wgs <- "data/WGS.hg38.strelka2.unpadded.bed"
opt$bed_wxs <- "data/WXS.hg38.100bp_padded.bed"
opt$annot_rds <- "analyses/snv-callers/ref_files/hg38_genomic_region_annotation.rds"
opt$no_region <- FALSE
opt$overwrite <- TRUE

# Coerce to numeric
opt$vaf_filter <- as.numeric(opt$vaf_filter)

# Make everything relative to root path
opt$maf <- file.path(root_dir, opt$maf)
opt$sql_file <- file.path(root_dir, opt$sql_file)
opt$bed_wgs <- file.path(root_dir, opt$bed_wgs)
opt$bed_wxs <- file.path(root_dir, opt$bed_wxs)

# Bring along the file suffix. Make to lower.
file_suffix <- tolower(opt$file_save)

# Check that the file format is supported
if (!(file_suffix %in% c("rds", "tsv", "none"))) {
  warning("Option used for file format (-f) is not supported. Only 'tsv' or 'rds'
          files are supported. Defaulting to not save any file besides the sqlite
          file.")
  opt$file_save <- "none"
  file_suffix <- "none"
}

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$maf, opt$sql_file, opt$bed_wgs, opt$bed_wxs
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
# Specify file names if options are set to save to rds or tsv files
if (opt$file_save != "none") {
  vaf_file <- file.path(caller_results_dir, paste0(opt$label, "_vaf.", file_suffix))
  region_annot_file <- file.path(caller_results_dir, paste0(opt$label, "_region.", file_suffix))
  tmb_file <- file.path(caller_results_dir, paste0(opt$label, "_tmb.", file_suffix))
  
  # Declare metadata file name for this caller
  metadata_file <- file.path(
    caller_results_dir,
    paste0(opt$label, "_metadata_filtered.", file_suffix)
  )
}

################## Calculate VAF and set up other variables ####################
# Print out warning if this file is going to be overwritten
if (file.exists(vaf_file)) {
  warning("Overwriting existing VAF file.")
}
# Print out progress message
message(paste("Calculating, sampling, and merging VAF for", opt$label, "MAF data..."))

# Print progress message
message(paste("Reading in", opt$maf, "MAF data..."))

# Use the premade function to calculate VAF this will also merge the metadata
vaf_df <- set_up_maf(opt$maf, opt$sql_file, opt$label, opt$overwrite,
                     vaf_cutoff = opt$vaf_filter)

message(paste(nrow(vaf_df), "mutations left after filter and merge"))

# Only if the file is being saved besides the SQL files
if (opt$file_save != "none") {
  # If the file exists or the overwrite option is not being used, calculate VAF
  if (file.exists(vaf_file) && !opt$overwrite) {
    # Stop if this file exists and overwrite is set to FALSE
    warning(cat(
      "The VAF file already exists: \n",
      vaf_file, "\n",
      "Use --overwrite if you want to overwrite it."
    ))
  } else {
    # Write the vaf file
    if (opt$file_save == "tsv") {
      vaf_df %>%
        readr::write_tsv(vaf_file)
    }
    if (opt$file_save == "rds") {
      vaf_df %>%
        readr::write_rds(vaf_file)
    }
    # Print out saved message
    message(paste("VAF calculations saved to: \n", vaf_file))
  }
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
    maf_annot <- annotr_maf(vaf_df, annotation_file = opt$annot_rds)

 if (opt$file_save != "none") {
    # Write the region file
    if (opt$file_save == "tsv") {
      maf_annot %>%
        readr::write_tsv(region_annot_file)
    }
    if (opt$file_save == "rds") {
      maf_annot %>%
        readr::write_rds(region_annot_file)
    }
      # Print out completion message
      message(paste("Region annotation file saved to:", region_annot_file))
    }
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
  if (opt$file_save == "tsv") {
    tmb_df %>%
      readr::write_tsv(tmb_file)
  }
  if (opt$file_save == "rds") {
    tmb_df %>%
      readr::write_rds(tmb_file)
  }
  
  if (opt$file_save != "none") {
    # Print out completion message
    message(paste("TMB calculations saved to:", tmb_file))
  }
}
