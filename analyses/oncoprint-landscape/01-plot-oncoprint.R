# This script displays an oncoprint displaying the landscape across PBTA given
# the relevant metadata and output MAF files from the snv callers strelka2,
# mutect2, lancet, and vardict. It addresses issue #6 in the OpenPBTA-analysis
# github repository.
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

# Define a color vector for plots
colores <- c(
  "Missense_Mutation" = "#35978f",
  "Nonsense_Mutation" = "#000000",
  "Frame_Shift_Del" = "#56B4E9",
  "Frame_Shift_Ins" = "#FFBBFF",
  "Splice_Site" = "#F0E442",
  "Translation_Start_Site" = "#191970",
  "Nonstop_Mutation" = "#545454",
  "In_Frame_Del" = "#CAE1FF",
  "In_Frame_Ins" = "#FFE4E1",
  "Stop_Codon_Ins" = "#CC79A7",
  "Start_Codon_Del" = "#56B4E9",
  "Fusion" = "#7B68EE",
  "Multi_Hit" = "#f46d43",
  "Hom_Deletion" = "#313695",
  "Hem_Deletion" = "#abd9e9",
  "Amplification" = "#c51b7d",
  "Loss" = "#0072B2",
  "Gain" = "#D55E00",
  "High_Level_Gain" = "#FF0000",
  "Multi_Hit_Fusion" = "#CD96CD"
)

#### Functions -----------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

plot_oncoplot <- function(maf_object, filename) {
  # Given a maf file and a filename, plot an oncoprint of the variants in the
  # dataset, save as a png file and display.
  # Args:
  #   maf_object: name or path to a maf object
  #   filename: name to save the png file as
  # Return:
  #   oncoprint: plot produced using the maftools `oncoplot` function
  
  # Plot and save the oncoprint
  png(
    file.path(plots_dir, filename),
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
    colors = colores
  )
  dev.off()
}

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

genes <- NULL
cnv_file <- NULL

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
    action = "store_true",
    type = "character",
    default = NULL,
    help = "file path to SEG file that contains cnv information"
  ),
  optparse::make_option(
    c("-f", "--fusion_file"),
    action = "store_true",
    type = "character",
    default = NULL,
    help = "file path to file that contains fusion information"
  ),
  optparse::make_option(
    c("-g", "--goi_list"),
    action = "store_true",
    type = "character",
    default = NULL,
    help = "file path to file that contains list of genes to include on
            oncoprint"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

maf <- opt$maf_file

if (!is.null(opt$cnv_file)) {
  cnv <- opt$cnv_file
}

if (!is.null(opt$fusion_file)) {
  fusion <- opt$fusion_file
}

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
  cnv_file <- data.table::fread(cnv, stringsAsFactors = FALSE)
}
# Read in fusion file
if (!is.null(opt$fusion_file)) {
  fusion_file <- data.table::fread(fusion, stringsAsFactors = FALSE)
}

#### Incorporate Fusion Data ---------------------------------------------------

if (!is.null(opt$fusion_file)) {
  fus <- fusion_file
  
  fus <- fus %>%
    dplyr::rename(Fusion = "gene1--gene2")
  
  # Separate fusion gene partners and add variant classification and center
  fus_sep <- fus %>%
    tidyr::separate(Fusion, c("5'-gene", "3'-gene"), sep = "--")
  
  # Determine which samples have mult-fused fusions
  multi <- fus_sep[, c("tumor_id", "5'-gene", "3'-gene")]
  multi$ID <- paste0(multi$tumor_id, ";", multi$`5'-gene`)
  reformat_multi <- multi %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(tumor_id = dplyr::n()) %>%
    as.data.frame()
  
  reformat_multi$Variant_Classification <-
    ifelse(reformat_multi$tumor_id == 1, "Fusion", "Multi_Hit_Fusion")
  reformat_multi <- reformat_multi %>%
    tidyr::separate(ID, c("Tumor_Sample_Barcode", "Hugo_Symbol"), sep = ";")
  reformat_multi$tumor_id <- NULL
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
    removeDuplicatedVariants = F,
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
    read.delim(opt$goi_list,
               sep = "\t",
               header = F,
               as.is = T)
  
  # Get top mutated this data and goi list
  gene_sum <- mafSummary(maf_object)$gene.summary
  
  # Subset for genes in the histology-specific list
  subset_gene_sum <- subset(gene_sum, Hugo_Symbol %in% goi_list$V1)
  
  # Get top altered genes
  goi_ordered <-
    subset_gene_sum[order(subset_gene_sum$AlteredSamples, decreasing = T),]
  
  # Select N top genes
  N <- ifelse(nrow(goi_ordered) > 20, 20, nrow(goi_ordered))
  goi_ordered_N <- goi_ordered[1:N,]
  genes <- goi_ordered_N$Hugo_Symbol
}

#### Plot and Save Oncoprint ---------------------------------------------------

plot_oncoplot(maf_object, "maf_oncoprint.png")
