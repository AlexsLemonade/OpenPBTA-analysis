# Run variant caller evaluation for a given MAF file.
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# -label : Label to be used for folder and all output. eg. 'strelka2'. Default is 'maf'.
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

# Coerce to numeric
opt$vaf_filter <- as.numeric(opt$vaf_filter)

# Make everything relative to root path
opt$maf <- file.path(root_dir, opt$maf)
opt$sql_file <- file.path(root_dir, opt$sql_file)
opt$bed_wgs <- file.path(root_dir, opt$bed_wgs)
opt$bed_wxs <- file.path(root_dir, opt$bed_wxs)

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

################## Calculate VAF and set up other variables ####################
# Print out progress message
message(paste("Calculating, sampling, and merging VAF for", opt$label, "MAF data..."))
  
# Use the premade function to calculate VAF this will also merge the metadata
set_up_maf(opt$maf, opt$sql_file, opt$label, opt$overwrite, 
           vaf_cutoff = opt$vaf_filter)

# Print out saved message
message(paste("VAF calculations saved to: \n", opt$sql_file))

######################### Annotate genomic regions #############################
if (opt$no_region) {
  # Print out progress message
  message(paste("Annotating genomic regions for", opt$label, "MAF data..."))

  # Annotation genomic regions
  annotr_maf(annotation_file = opt$annot_rds)
  
  # Print out completion message
  message(paste("Region annotation file saved to:", opt$sql_file))
}

############################# Calculate TMB ####################################
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

# Filter out mutations for WXS that are outside of these BED regions.
wxs_bed_filter(wxs_bed_file = opt$bed_wxs)

# Calculate TMBs and write to SQL file
calculate_tmb(
  wgs_size = wgs_genome_size,
  wxs_size = wxs_exome_size
)

# Print out completion message
message(paste("TMB calculations saved to:", opt$sql_file))
