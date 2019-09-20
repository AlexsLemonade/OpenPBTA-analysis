# Run variant caller evaluation for a given MAF file.
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# -m : Path to MAF file to be analyzed. Can be .gz compressed.
# -b : File path that specifies the caller specific WGS BED regions file that is
#      saved as a TSV but has CHR, start, and stop as it's first three columns
#      Note that column names are ignored.
# -l : Label to be used for folder and all output. eg. 'strelka2'. Optional.
#      Default is 'maf'
#
# Command line example:
#
# Rscript 01-calculate_vaf_tmb.R \
# -m ../../data/pbta-snv-strelka2.vep.maf.gz \
# -b bed_regions/strelka2_wgs_bed_regions.tsv \
# -l strelka2
#
# Run set up script
source(file.path("scripts", "00-set_up.R"))

# Load library:
library(optparse)

################################ Set up options ################################

# Set up optparse options
option_list <- list(
  make_option(opt_str = c("-m", "--maf"), type = "character", default = "none",
              help = "Path to MAF file to be analyzed. Can be .gz compressed.",
              metavar = "character"),
  make_option(opt_str = c("-b", "--bed_wgs"), type = "character",
              default = "none", help = "File path that specifies the caller specific
                BED regions file.",
              metavar = "character"),
  make_option(opt_str = c("-l", "--label"), type = "character",
              default = "maf", help = "Label to be used for folder and all
                output. eg. 'strelka2'. Optional. Default is 'maf'",
              metavar = "character")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Stop if no input data is specified
if (opt$maf == "none") {
  stop("Error: no specified MAF file to analyze at this file path. Use -m
         to give a file.path")
}

# Stop if no input data is specified
if (opt$bed_wgs == "none") {
  stop("Error: no specified BED file to analyze at this file path. Use -b
         to give a file.path")
}

################## Create output directories for this caller ###################
caller_results_dir <- file.path(base_results_dir, opt$label)

# Make caller specific results folder
if (!dir.exists(caller_results_dir)) {
  dir.create(caller_results_dir)
}

####################### File paths for files we will create ####################
vaf_file <- file.path(caller_results_dir, paste0(opt$label, "_vaf.tsv"))
region_annot_file <- file.path(caller_results_dir, paste0(opt$label, "_region.tsv"))
tmb_file <- file.path(caller_results_dir, paste0(opt$label, "_tmb.tsv"))

# Declare metadata file name for this caller
metadata_file <- file.path(caller_results_dir, 
                           paste0(opt$label, "_metadata_filtered.tsv"))

########################### Set up this caller's data ##########################
# Print progress message
message(paste("Reading in", opt$maf, "MAF data..."))

# Read in this MAF, skip the version number
# maf_df <- data.table::fread(opt$maf, skip = 1, data.table = FALSE)
maf_df <- data.table::fread(opt$maf, data.table = FALSE)

# Print progress message
message(paste("Setting up", opt$label, "metadata..."))

# Isolate metadata to only the samples that are in the datasets
metadata <- readr::read_tsv(original_metadata) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% maf_df$Tumor_Sample_Barcode) %>%
  dplyr::distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE) %>%
  dplyr::arrange() %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) %>%
  readr::write_tsv(metadata_file)

# Make sure that we have metadata for all these samples.
if (!all(unique(maf_df$Tumor_Sample_Barcode) %in% metadata$Tumor_Sample_Barcode)) {
  stop("There are samples in this MAF file that are not in the metadata.")
}

# Tack on the metadata so we have this info
maf_df <- maf_df %>%
  dplyr::left_join(metadata, by = "Tumor_Sample_Barcode") %>%
  # Get rid of any variables that have completely NAs.
  dplyr::select(-which(apply(is.na(.), 2, all)))

################## Calculate VAF and set up other variables ####################
message(paste("Calculating VAF for", opt$label, "MAF data..."))

# Use the premade function to calculate VAF
vaf_df <- calculate_vaf(maf_df) %>%
  readr::write_tsv(vaf_file)

######################### Annotate genomic regions #############################
message(paste("Annotating genomic regions for", opt$label, "MAF data..."))

# Annotation genomic regions
maf_annot <- annotr_maf(vaf_df) %>%
  readr::write_tsv(region_annot_file)

############################# Calculate TMB ####################################
message(paste("Calculating TMB for", opt$label, "MAF data..."))

### Set up BED region files for TMB calculations
wgs_bed <- readr::read_tsv(opt$bed_wgs, col_names = FALSE)
wxs_bed <- readr::read_tsv(wxs_bed_file, col_names = FALSE)

# Calculate size of genome surveyed
wgs_genome_size <- sum(wgs_bed[, 3] - wgs_bed[, 2])
wxs_exome_size <- sum(wxs_bed[, 3] - wxs_bed[, 2])

# Print out these genome sizes
cat(
  "WGS size in bp:", wgs_genome_size,
  "WXS size in bp:", wxs_exome_size
)

# Filter out mutations for WXS that are outside of these BED regions.
maf_wxs_filtered <- wxs_bed_filter(vaf_df, wxs_bed)

# Calculate TMBs and write to TMB file
tmb_df <- calculate_tmb(maf_wxs_filtered,
  wgs_size = wgs_genome_size,
  wxs_size = wxs_exome_size
) %>%
  readr::write_tsv(tmb_file)
