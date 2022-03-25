# This script outputs the oncoprint N counts table 

#### Set Up --------------------------------------------------------------------

# Load libraries
library(dplyr)
library(maftools)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to output directory for table produced
results_dir <-
  file.path(root_dir, "analyses", "oncoprint-landscape", "tables")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

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
    c("-o", "--output_file"),
    type = "character",
    default = NULL,
    help = "output filename, not the full path"
  ),
  optparse::make_option(
    c("--include_introns"),
    action = "store_true",
    default = FALSE,
    help = "logical statement on whether to include intronic variants in oncoprint plot"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

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

# Read in fusion file and join
fusion_df <- readr::read_tsv(opt$fusion_file)

#### Set up oncoprint counts ---------------------------------------------------

# Broad histologies that fall under "Other CNS" for oncoprint purposes
other_cns_histologies <- c(
  "Ependymal tumor",
  "Tumors of sellar region",
  "Neuronal and mixed neuronal-glial tumor",
  "Tumor of cranial and paraspinal nerves",
  "Meningioma",
  "Mesenchymal non-meningothelial tumor",
  "Germ cell tumor",
  "Choroid plexus tumor",
  "Histiocytic tumor",
  "Tumor of pineal region",
  "Metastatic tumors",
  "Other astrocytic tumor",
  "Lymphoma",
  "Melanocytic tumor",
  "Other tumor"
)

# All broad histologies we want to generate sample counts for
broad_histologies_included <- c(
  other_cns_histologies,
  "Low-grade astrocytic tumor",
  "Embryonal tumor",
  "Diffuse astrocytic and oligodendroglial tumor"
)

# Sample IDs to include
tumor_barcodes_to_include <- metadata %>%
  filter(broad_histology %in% broad_histologies_included) %>%
  pull(Tumor_Sample_Barcode)

# Filter data files to only include the relevant samples & then pull the IDs
# We're only counting IDs that are present in the relevant file using this
# approach.
maf_filtered_ids <- maf_df %>%
  filter(Tumor_Sample_Barcode %in% tumor_barcodes_to_include) %>%
  pull(Tumor_Sample_Barcode)

# CNV alterations
cnv_filtered_ids <- cnv_df %>%
  filter(Tumor_Sample_Barcode %in% tumor_barcodes_to_include) %>%
  pull(Tumor_Sample_Barcode)

# Fusion alterations
fusion_filtered_ids <- fusion_df %>%
  filter(Tumor_Sample_Barcode %in% tumor_barcodes_to_include) %>%
  pull(Tumor_Sample_Barcode)

# Take the union of filtered IDs
filtered_ids <- union(maf_filtered_ids, cnv_filtered_ids ) %>% 
  union(fusion_filtered_ids)

# Create the table of counts
sample_counts_df <- metadata %>%
  filter(Tumor_Sample_Barcode %in% filtered_ids,
         !is.na(broad_histology)) %>%
  select(Tumor_Sample_Barcode, broad_histology) %>%
  mutate(broad_histology = case_when(
    broad_histology %in% other_cns_histologies ~ "Other CNS",
    TRUE ~ broad_histology
  )) %>%
  distinct() %>%
  count(broad_histology)

# Write to file
output_file_path <- file.path(results_dir, opt$output_file)
readr::write_tsv(sample_counts_df, output_file_path)
