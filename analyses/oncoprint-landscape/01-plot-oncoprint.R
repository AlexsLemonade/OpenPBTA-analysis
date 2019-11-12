# This script displays an oncoprint displaying the landscape across PBTA given
# the relevant metadata. It addresses issue #6 in the OpenPBTA-analysis
# github repository. It uses the output of 00-map-to-sample_id.R.
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
    help = "file path to file that contains list of genes to include on
            oncoprint"
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

# Define cnv_file object here as it still needs to be defined for the `read.maf`
# function, even if it is NULL
cnv_file <- opt$cnv_file

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
  cnv_file <- readr::read_tsv(cnv_file)
}

# Read in fusion file and join
if (!is.null(opt$fusion_file)) {
  fusion_file <- readr::read_tsv(opt$fusion_file)
  maf_df <- dplyr::bind_rows(maf_df, fusion_file)
}

#### Convert into MAF object ---------------------------------------------------

maf_object <-
  read.maf(
    maf = maf_df,
    clinicalData = metadata,
    cnTable = cnv_file,
    removeDuplicatedVariants = FALSE,
    vc_nonSyn = c(
      "Frame_Shift_Del",
      "Frame_Shift_Ins",
      "Splice_Site",
      "Nonsense_Mutation",
      "Nonstop_Mutation",
      "In_Frame_Del",
      "In_Frame_Ins",
      "Missense_Mutation",
      "Fusion",
      "Multi_Hit",
      "Multi_Hit_Fusion",
      "Hom_Deletion",
      "Hem_Deletion",
      "amplification",
      "gain",
      "loss"
    )
  )

#### Specify genes -------------------------------------------------------------

if (!is.null(opt$goi_list)) {
  # Read in gene list
  goi_list <-
    read.delim(
      file.path(opt$goi_list),
      sep = "\t",
      header = FALSE,
      as.is = TRUE
    )

  # Get top mutated this data and goi list
  gene_sum <- mafSummary(maf_object)$gene.summary

  # Subset for genes in the histology-specific list
  subset_gene_sum <- subset(gene_sum, Hugo_Symbol %in% goi_list$V1)

  # Get top altered genes
  goi_ordered <-
    subset_gene_sum[order(subset_gene_sum$AlteredSamples, decreasing = TRUE), ]

  # Select n top genes
  num_genes <- ifelse(nrow(goi_ordered) > 20, 20, nrow(goi_ordered))
  goi_ordered_num <- goi_ordered[1:num_genes, ]
  genes <- goi_ordered_num$Hugo_Symbol

} else {
  # If a gene list is not supplied, we do not want the `oncoplot` function to
  # filter the genes to be plotted, so we assign NULL to the `genes` object.
  genes <- NULL
}

#### Plot and Save Oncoprint ---------------------------------------------------

# Given a maf file, plot an oncoprint of the variants in the
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
  genes = genes,
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
