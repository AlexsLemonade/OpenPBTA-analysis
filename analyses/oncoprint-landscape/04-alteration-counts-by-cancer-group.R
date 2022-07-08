# J. Taroni for ALSF CCDL 2022
# Counts alterations per cancer group to be reported in the manuscript text.
# These are counts for oncoprint plots that use genes of interest lists.
# The mappings between cancer groups and genes of interest lists in available in
# data/cancer_group_goi_list_mapping.tsv in this module, which is read in 
# directly, rather than passed in via an argument. Data files need to have been
# processed with 01-map-to-sample-id.R. The resulting tables are output into
# tables/cancer_group_counts/<subdirectory>, where <subdirectory> is set via an
# argument.
#
# USAGE: 
#  Rscript --vanilla 04-alteration-counts-by-cancer-group.R \
#    --maf_file ../../scratch/oncoprint_files/primary_only_maf.tsv \
#    --cnv_file ../../scratch/oncoprint_files/primary_only_cnv.tsv \
#    --fusion_file ../../scratch/oncoprint_files/primary_only_fusions.tsv \
#    --metadata_file ../../data/pbta-histologies.tsv \
#    --palette_file ../../figures/palettes/broad_histology_cancer_group_palette.tsv \
#    --subdirectory primary_only


#### Set Up --------------------------------------------------------------------

# Load libraries
library(dplyr)
library(maftools)
library(stringr)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# This module
module_dir <- file.path(root_dir, "analyses", "oncoprint-landscape")

# Source the custom functions script
source(
  file.path(
    module_dir,
    "util",
    "oncoplot-functions.R"
  )
)

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-m", "--maf_file"),
    type = "character",
    default = NULL,
    help = "file path to MAF file that contains SNV information",
  ),
  optparse::make_option(
    c("-c", "--cnv_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains CNV information"
  ),
  optparse::make_option(
    c("-f", "--fusion_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains fusion information"
  ),
  optparse::make_option(
    c("-s", "--metadata_file"),
    type = "character",
    default = NULL,
    help = "file path to the histologies file"
  ),
  optparse::make_option(
    c("-p", "--palette_file"),
    type = "character",
    default = NULL,
    help = "file path to the palette file"
  ),    
  optparse::make_option(
    c("-d", "--subdirectory"),
    type = "character",
    default = NULL,
    help = "subdirectory, used in practice to partition primary from primary-plus"
  ),
  optparse::make_option(
    c("--include_introns"),
    action = "store_true",
    default = FALSE,
    help = "logical statement on whether to include intronic variants in MAF file"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Path to output directory for tables produced
results_dir <- file.path(module_dir, 
                         "tables", 
                         "cancer_group_counts",
                         opt$subdirectory)

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <- readr::read_tsv(opt$metadata_file, guess_max = 10000) %>%
  rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF file
maf_df <- data.table::fread(opt$maf_file,
                            stringsAsFactors = FALSE,
                            data.table = FALSE)

if (!opt$include_introns) {
  maf_df <- maf_df %>%
    filter(Variant_Classification != "Intron")
}

# Read in cnv file
cnv_df <- readr::read_tsv(opt$cnv_file) 

# Read in fusion file
fusion_df <- readr::read_tsv(opt$fusion_file)

# Read in palette file, and join with the metadata file
palette_df <- readr::read_tsv(opt$palette_file)
metadata <- left_join(
  select(palette_df, cancer_group, cancer_group_display, broad_histology), 
  metadata
)

#### Read in cancer group to goi list mapping ----------------------------------

goi_file <- file.path(module_dir,
                      "data",
                      "cancer_group_goi_list_mapping.tsv")
cancer_group_goi_df <- readr::read_tsv(goi_file)

for (cg in cancer_group_goi_df$cancer_group_display) {
  
  # Filtered metadata by cancer group
  filtered_metadata_df <- metadata %>%
    filter(cancer_group_display == cg)
  
  # Get list of samples
  cg_samples <- filtered_metadata_df %>%
    pull(Tumor_Sample_Barcode)
  
  # Filter all the data files to only samples in this cancer group
  filtered_maf_df <- maf_df %>%
    filter(Tumor_Sample_Barcode %in% cg_samples)
  filtered_cnv_df <- cnv_df %>%
    filter(Tumor_Sample_Barcode %in% cg_samples)
  filtered_fusion_df <- fusion_df %>%
    filter(Tumor_Sample_Barcode %in% cg_samples)
  
  # We need to read in the GOI to get the correct counts
  cg_goi_path <- cancer_group_goi_df %>%
    filter(cancer_group_display == cg) %>%
    # Have to ensure only 1 file is returned
    select(goi_list_file) %>%
    distinct() %>%
    pull(goi_list_file)
  
  goi_list <- readr::read_tsv(file.path(module_dir, cg_goi_path)) %>%
    as.matrix()
  
  # Sometimes we get
  # Error in read.maf(maf = maf_df, clinicalData = metadata, cnTable = cnv_df,  :
  # No non-synonymous mutations found
  # Thanks SO: https://stackoverflow.com/questions/49966585/using-trycatch-function-of-r-in-a-loop
  cg_maf_object <- tryCatch(
    prepare_maf_object(
      maf_df = filtered_maf_df,
      cnv_df = filtered_cnv_df,
      metadata = filtered_metadata_df,
      fusion_df = filtered_fusion_df
    ),
    error = function(e) { return(0) }
  )
  
  # If prepare_maf_object() above threw an error, cg_maf_object will be 0 and
  # we want to skip to the next cancer group
  if (is.numeric(cg_maf_object)) next 

  # Sometimes there will be no non-synonymous variants and subsetMaf() will
  # throw an error
  filtered_maf_object <- tryCatch(
    subsetMaf(
      maf = cg_maf_object,
      tsb = filtered_metadata_df$Tumor_Sample_Barcode,
      genes = goi_list,
      mafObj = TRUE
    ),
    error = function(e) { return(0) }
  )

  # If subsetMaf() above threw an error, filtered_maf_object will be 0 and
  # we want to skip to the next cancer group
  if (is.numeric(filtered_maf_object)) next 

  # Get top mutated genes per this subset object
  gene_sum <- mafSummary(filtered_maf_object)$gene.summary

  # Write this to file
  # First we'd like a "nice" version of cancer group for the file name
  cg_nice <- str_replace_all(  # Make spaces dashes
    str_replace_all(  # Remove all apostrophe
      str_to_lower(cg),  # Make all lower case
      "'", ""), 
    " ", "-")

  # Now construct the path and write to file
  output_file <- file.path(results_dir,
                           paste0(cg_nice, "_oncoprint_alteration_counts.tsv"))
  readr::write_tsv(gene_sum, output_file)

}
