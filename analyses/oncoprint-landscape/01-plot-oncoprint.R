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
    c("-b", "--broad_histology"),
    type = "character",
    default = NULL,
    help = "optional name of `broad_histology` value to plot associated oncoprint"
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
metadata <- readr::read_tsv(opt$metadata_file, guess_max = 10000) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF file
maf_df <- data.table::fread(opt$maf_file,
                            stringsAsFactors = FALSE,
                            data.table = FALSE)

# Read in cnv file
if (!is.null(opt$cnv_file)) {
  cnv_df <- readr::read_tsv(opt$cnv_file) %>%
    dplyr::mutate(Variant_Classification = dplyr::case_when(Variant_Classification == "loss" ~ "Del",
                                                            Variant_Classification %in% c("gain", "amplification") ~ "Amp",
                                                            TRUE ~ as.character(Variant_Classification)))
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

#### Set up oncoprint annotation objects --------------------------------------
# Read in histology standard color palette for project
histology_label_mapping <- readr::read_tsv(
  file.path(root_dir,
            "figures",
            "palettes", 
            "histology_label_color_table.tsv")) %>%
  # Select just the columns we will need for plotting
  dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes)

# Join on these columns to the metadata
metadata <- metadata %>% 
  dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID") %>% 
  # Reorder display_group based on display_order
  dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))

# Filter to the metadata associated with the broad histology value, if provided
if (!is.null(opt$broad_histology)) {
  
  if (!opt$broad_histology == "Other CNS") {
    metadata <- metadata %>%
      dplyr::filter(broad_histology == opt$broad_histology)
  } else {
    metadata <- metadata %>%
      dplyr::filter(
        broad_histology %in% c(
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
      )
  }
  
  # Now filter the remaining data files
  maf_df <- maf_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)
  
  cnv_df <- cnv_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)
  
  fusion_df <- fusion_df %>%
    dplyr::filter(Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode)
  
}

# Read in the oncoprint color palette
oncoprint_col_palette <- readr::read_tsv(file.path(
  root_dir,
  "figures",
  "palettes",
  "oncoprint_color_palette.tsv"
)) %>%
  # Use deframe so we can use it as a recoding list
  tibble::deframe()

# Color coding for `display_group` classification
# Get unique tumor descriptor categories
histologies_color_key_df <- metadata %>%
  dplyr::arrange(display_order) %>%
  dplyr::select(display_group, hex_codes) %>%
  dplyr::distinct()

# Make color key specific to these samples
histologies_color_key <- histologies_color_key_df$hex_codes
names(histologies_color_key) <- histologies_color_key_df$display_group

# Now format the color key objet into a list
annotation_colors <- list(display_group = histologies_color_key)

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
  file.path(plots_dir, tolower(gsub(" ", "-", opt$png_name))),
  width = 65,
  height = 30,
  units = "cm",
  res = 300
)
oncoplot(
  maf_object,
  clinicalFeatures = "display_group",
  genes = goi_list,
  logColBar = TRUE,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = TRUE,
  removeNonMutated = TRUE,
  annotationFontSize = 1.0,
  SampleNamefontSize = 0.5,
  fontSize = 0.7,
  colors = oncoprint_col_palette,
  annotationColor = annotation_colors,
  bgCol = "#F5F5F5",
  top = 25
)
dev.off()
