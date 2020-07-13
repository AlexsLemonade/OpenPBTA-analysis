# This script displays an oncoprint displaying the landscape across PBTA given
# the relevant metadata. It addresses issue #6 in the OpenPBTA-analysis
# github repository. It uses the output of 00-map-to-sample_id.R. It can
# accept a gene list file or a comma-separated set of gene list files that will
# be concatenated and restricts plotting to those genes (via --goi_list).
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

# Load maftools
library(maftools)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to output directory for plots produced
plots_dir <-
  file.path(root_dir, "analyses", "oncoprint-landscape", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

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
    c("-g", "--goi_list"),
    type = "character",
    default = NULL,
    help = "comma-separated list of genes of interest files that contain the
            genes to include on oncoprint"
  ),
  optparse::make_option(
    c("-p", "--png_name"),
    type = "character",
    default = NULL,
    help = "oncoprint output png file name"
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
goi_list <- opt$goi_list

#### Functions ----------------------------------------------------------------

read_genes <- function(gene_list) {
  # This function takes in the file path to a gene list and pulls out
  # the gene information from that list
  #
  # Args:
  #   gene_list: file path to genes of interest file
  #
  # Return:
  #   genes: a vector of genes from the genes of interest file

  genes <- readr::read_tsv(gene_list) %>%
    dplyr::pull("gene")
}

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <- readr::read_tsv(opt$metadata_file,
                            guess_max = 10000) %>%
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

# Read in gene information from the list of genes of interest files
if (!is.null(opt$goi_list)) {
  goi_files <- unlist(stringr::str_split(goi_list, ",| "))
  # Read in using the `read_genes` custom function and unlist the gene column
  # data from the genes of interest file paths given
  goi_list <- lapply(goi_files, read_genes)
    # Include only the unique genes of interest
  goi_list <- unique(unlist(goi_list))
}

# Read in recurrent focal CNVs file
if (!is.null(opt$focal_file)) {
  focal_df <- readr::read_tsv(file.path(opt$focal_file))
}

# Read in histology standard color palette for project
histology_col_palette <-
  readr::read_tsv(file.path(
    root_dir,
    "figures",
    "palettes",
    "histology_color_palette.tsv"
  ))

# Read in the oncoprint color palette
oncoprint_col_palette <- readr::read_tsv(file.path(
  root_dir,
  "figures",
  "palettes",
  "oncoprint_color_palette.tsv"
)) %>%
  # Use deframe so we can use it as a recoding list
  tibble::deframe()

#### Set up oncoprint annotation objects --------------------------------------

# Color coding for `short_histology` classification
# Get unique tumor descriptor categories
short_histologies <- unique(metadata$short_histology) %>%
  tidyr::replace_na("none") %>%
  sort()

# Save the vector of hex codes from the short histology palette
short_histology_col_key <- histology_col_palette$hex_codes

# Now assign the color names
names(short_histology_col_key) <- short_histologies

# Now format the color key objet into a list
annotation_colors <- list(short_histology = short_histology_col_key)

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
  clinicalFeatures = "short_histology",
  genes = goi_list,
  logColBar = TRUE,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = TRUE,
  removeNonMutated = TRUE,
  annotationFontSize = 1.0,
  SampleNamefontSize = 0.5,
  fontSize = 0.7,
  colors = oncoprint_col_palette,
  annotationColor = annotation_colors
)
dev.off()
