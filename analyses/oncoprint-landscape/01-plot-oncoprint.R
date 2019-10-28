# This script displays an oncoprint displaying the landscape across PBTA given
# the relevant metadata and output MAF files from the snv callers strelka2,
# mutect2, lancet, and vardict. It addresses issue #6 in the OpenPBTA-analysis
# github repository.
#
# Code adapted from the PPTC PDX Oncoprint Generation repository here:
# https://github.com/marislab/create-pptc-pdx-oncoprints/tree/master/R
#
# Chante Bethell for CCDL 2019 and Jo Lynne Rokita
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/oncoprint-landscape/01-plot-oncoprint.R'

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
    help = "file path to MAF file that contains snv information",
  ),
  optparse::make_option(
    c("-c", "--cnv_file"),
    type = "character",
    default = NULL,
    help = "file path to SEG file that contains cnv information"
  ),
  optparse::make_option(
    c("-f", "--fusion_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains fusion information"
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

maf <- opt$maf_file

# Define cnv_file object here as it still needs to be defined for the `read.maf`
# function, even if it is NULL
cnv_file <- opt$cnv_file

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <-
  readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"))

# Rename for maftools function
metadata <- metadata %>%
  dplyr::rename("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")

# Read in MAF file
maf_df <- data.table::fread(maf, stringsAsFactors = FALSE)

# Read in cnv file
if (!is.null(opt$cnv_file)) {
  cnv_file <- data.table::fread(opt$cnv_file, stringsAsFactors = FALSE)
  # TODO: Filter and set up `cnv_file` to be in the column format -
  # "Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification" as required by
  # the `read.maf function`
}

# Read in fusion file
if (!is.null(opt$fusion_file)) {
  fusion_file <- data.table::fread(opt$fusion_file, stringsAsFactors = FALSE)

  #### Incorporate Fusion Data -------------------------------------------------

  # Separate fusion gene partners and add variant classification and center
  fus_sep <- fusion_file %>%
    tidyr::separate(FusionName, c("5'-gene", "3'-gene"), sep = "--")

  # Determine which samples have multi-fused fusions
  multi <- fus_sep[, c("Sample", "5'-gene", "3'-gene")]
  multi$ID <- paste0(multi$Sample, ";", multi$`5'-gene`)
  reformat_multi <- multi %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(Sample = dplyr::n()) %>%
    as.data.frame()

  reformat_multi$Variant_Classification <-
    ifelse(reformat_multi$Sample == 1, "Fusion", "Multi_Hit_Fusion")
  reformat_multi <- reformat_multi %>%
    tidyr::separate(ID, c("Tumor_Sample_Barcode", "Hugo_Symbol"), sep = ";")
  reformat_multi$Sample <- NULL
  reformat_multi$Variant_Type <- "OTHER"

  # Merge with MAF
  maf_df <- dplyr::bind_rows(maf_df, reformat_multi)
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
      "Multi_Hit_Fusion"
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
  # filter the genes to be plotted, so we assign NULL to the `genes`` object.
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
  removeNonMutated = FALSE,
  annotationFontSize = 0.7,
  SampleNamefontSize = 0.5,
  fontSize = 0.7,
  colors = color_palette
)
dev.off()
