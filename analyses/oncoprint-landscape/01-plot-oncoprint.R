# This script displays an oncoprint displaying the landscape across PBTA given
# the relevant metadata. It addresses issue #6 in the OpenPBTA-analysis
# github repository. It uses the output of 00-map-to-sample_id.R. It can
# accept a gene list (via --goi_list) that restricts plotting to those genes.
# Gene lists should be a single column with no header and use gene symbols.
#
# Code adapted from the PPTC PDX Oncoprint Generation repository here:
# https://github.com/marislab/create-pptc-pdx-oncoprints/tree/master/R
#
# Chante Bethell for CCDL 2019 and Jo Lynne Rokita
#
# EXAMPLE USAGE:
#
# Rscript --vanilla 01-plot-oncoprint.R \
#  --maf_file ../../scratch/all_primary_samples_maf.tsv \
#  --cnv_file ../../scratch/all_primary_samples_cnv.tsv \
#  --fusion_file ../../scratch/all_primary_samples_fusions.tsv \
#  --metadata_file ../../data/pbta-histologies.tsv \
#  --goi_list ${genes_list} \
#  --png_name ${primary_filename}_goi_oncoprint.png


#### Set Up --------------------------------------------------------------------

# Install maftools
if (!("maftools" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("maftools")
}
library(maftools)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to the data obtained via `bash download-data.sh`.
data_dir <- file.path(root_dir, "data")

# Path to output directory for plots produced
plots_dir <-
  file.path(root_dir, "analyses", "oncoprint-landscape", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Source the color palette for plots
source(
  file.path(
    root_dir,
    "analyses",
    "oncoprint-landscape",
    "util",
    "oncoplot-palette.R"
  )
)

# Source the custom functions script
source(
  file.path(
    root_dir,
    "analyses",
    "oncoprint-landscape",
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
    c("-g", "--goi_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains list of genes to include on
            oncoprint"
  ),
  optparse::make_option(
    c("-p", "--png_name"),
    type = "character",
    default = NULL,
    help = "oncoprint output png file name"
  ),
  optparse::make_option(
    c("-o", "--focal_file"),
    type = "character",
    default = NULL,
    help = "file path to most focal CN units file"
  ),
  optparse::make_option(
    c("-r", "--recurrent_file"),
    type = "character",
    default = NULL,
    help = "file path to recurrent focal CN units file"
  )    
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Define cnv_file, fusion_file, and genes object here as they still need to
# be defined for the `prepare_and_plot_oncoprint` custom function (the 
# cnv_file specifically for the `read.maf` function within the custom function),
# even if they are NULL
cnv_df <- opt$cnv_file
fusion_df <- opt$fusion_file
goi_list <- opt$goi_file

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <- readr::read_tsv(opt$metadata_file) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF file
maf_df <- data.table::fread(opt$maf_file,
                            stringsAsFactors = FALSE,
                            data.table = FALSE)

# Read in cnv file
if (!is.null(opt$cnv_file)) {
  cnv_df <- readr::read_tsv(opt$cnv_file)
}

# Read in fusion file and join
if (!is.null(opt$fusion_file)) {
  fusion_df <- readr::read_tsv(opt$fusion_file)
}

# Read in genes list
if (!is.null(opt$goi_file)) {
  goi_list <- readr::read_tsv(file.path(opt$goi_file)) %>%
    dplyr::pull("gene")
    
  # Filter `goi_list` to include only the unique genes of interest
  goi_list <- unique(goi_list)
}

#opt$focal_file <- file.path(root_dir, "analyses", "focal-cn-file-preparation", "results", "consensus_seg_most_focal_cn_status.tsv.gz")
# Read in most focal CN units file
if (!is.null(opt$focal_file)) {
most_focal_df <- readr::read_tsv(file.path(opt$focal_file))
}

#opt$recurrent_file <- file.path(root_dir, "analyses", "focal-cn-file-preparation", "results", "consensus_seg_recurrent_focal_cn_units.tsv")
if (!is.null(opt$recurrent_file)) {
# Read in recurrent focal CN units file
recurrent_focal_df <- readr::read_tsv(file.path(opt$recurrent_file))
}

#### Format recurrent focal CN object -----------------------------------------

if (!is.null(opt$focal_file) & !is.null(opt$recurrent_file)) {
# Let's filter the most focal calls to those that are considered recurrent, then
# filter to the primary plus samples
cnv_df <- most_focal_df %>%
  dplyr::filter(region %in% recurrent_focal_df$region,
         status != "uncallable") %>%
  dplyr::left_join(metadata, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Hugo_Symbol = region, Tumor_Sample_Barcode, Variant_Classification = status) %>%
  dplyr::filter(Tumor_Sample_Barcode %in% cnv_df$Tumor_Sample_Barcode)
}

#### Prepare MAF object for plotting ------------------------------------------

maf_object <- prepare_maf_object(
  maf_df = maf_df,
  cnv_df = cnv_df,
  metadata = metadata,
  fusion_df = fusion_df
)

#### Plot and Save Oncoprint --------------------------------------------------

# Given a maf object, plot an oncoprint of the variants in the
# dataset and save as a png file.
png(
  file.path(plots_dir, opt$png_name),
  width = 65,
  height = 30,
  units = "cm",
  res = 300
)
oncoplot(
  maf_object,
  clinicalFeatures = c(
    "broad_histology",
    "short_histology",
    "reported_gender",
    "tumor_descriptor",
    "molecular_subtype"
  ),
  genes = goi_list,
  logColBar = TRUE,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = TRUE,
  removeNonMutated = TRUE,
  annotationFontSize = 0.7,
  SampleNamefontSize = 0.5,
  fontSize = 0.7,
  colors = color_palette
)
dev.off()
