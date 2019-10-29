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
  ),
  optparse::make_option(
    c("-l", "--low_segmean_cutoff"),
    type = "double",
    default = 0.5, # Determined using CNVkit documentation (https://cnvkit.readthedocs.io/en/stable/calling.html)
    help = "low segmean cutoff to determine amplification and deletion"
  ),
  optparse::make_option(
    c("-s", "--high_segmean_cutoff"),
    type = "double",
    default = 1, # Determined using CNVkit documentation (https://cnvkit.readthedocs.io/en/stable/calling.html)
    help = "high segmean cutoff to distinguish hemizygous deletion from
            homozygous deletion"
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
  cnv_file <-
    data.table::fread(opt$cnv_file,
                      data.table = FALSE,
                      stringsAsFactors = FALSE)
  # TODO: Adapt `cnv_file` once the results of the consensus calls are obtained,
  # to be in the column format -
  # "Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification" as required by
  # the `read.maf function`

  # Create a dataframe with `Tumor_Sample_Barcode` and gene information from maf
  # dataframe
  select_maf_df <- maf_df %>%
    dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol)

  # Add gene information to the cnv dataframe
  cnv_df <- cnv_file %>%
    dplyr::inner_join(select_maf_df, by = c("ID" = "Tumor_Sample_Barcode")) %>%
    dplyr::distinct()

  # Define cutoffs for `Variant_Classification` column
  low_segmean_cutoff <- opt$low_segmean_cutoff
  high_segmean_cutoff <- opt$high_segmean_cutoff

  # Create `Variant_Classification` column
  cnv_df <- cnv_df %>%
    dplyr::mutate(
      Variant_Classification = dplyr::case_when(
        seg.mean < -low_segmean_cutoff &
          seg.mean > -high_segmean_cutoff ~ "Hem_Deletion",
        seg.mean < -high_segmean_cutoff ~ "Hom_Deletion",
        seg.mean > low_segmean_cutoff ~ "Amplification"
      )
    )

  # Make the `Variant_Classification` column a factor vector
  cnv_df$Variant_Classification <- as.factor(cnv_df$Variant_Classification)

  # Select the columns specified as input format for `read.maf` function
  cnv_df <- cnv_df %>%
    dplyr::select(Hugo_Symbol, ID, Variant_Classification)

  # Rename columns
  colnames(cnv_df) <-
    c("Hugo_Symbol",
      "Tumor_Sample_Barcode",
      "Variant_Classification")

  # Remove NA
  cnv_file <- cnv_df %>%
    dplyr::filter(!is.na(Variant_Classification))
}

# Read in fusion file
if (!is.null(opt$fusion_file)) {
  fusion_file <-
    data.table::fread(opt$fusion_file,
                      data.table = FALSE,
                      stringsAsFactors = FALSE)

  #### Incorporate Fusion Data -------------------------------------------------
  # TODO: Once the consensus calls of the fusion data are obtained, this section
  # will need to be adapted to the format of the fusion input file. For example,
  # the way we separate the genes out of `FusionName` may need to be adapted.

  # Separate fusion gene partners and add variant classification and center
  fus_sep <- fusion_file %>%
    # Separate the 5' and 3' genes
    tidyr::separate(FusionName, c("Gene1", "Gene2"), sep = "--") %>%
    dplyr::select(Sample, Gene1, Gene2)

  reformat_fusion <- fus_sep %>%
    # Here we want to tally how many times the 5' gene shows up as a fusion hit
    # in a sample
    dplyr::group_by(Sample, Gene1) %>%
    dplyr::tally() %>%
    # If the sample-5' gene pair shows up more than once, call it a multi hit
    # fusion
    dplyr::mutate(Variant_Classification =
                    dplyr::if_else(n == 1, "Fusion", "Multi_Hit_Fusion"),
                  # Required column for joining with MAF
                  Variant_Type = "OTHER") %>%
    # Correct format for joining with MAF
    dplyr::rename(Tumor_Sample_Barcode = Sample, Hugo_Symbol = Gene1) %>%
    dplyr::select(-n)

  # Merge with MAF
  maf_df <- dplyr::bind_rows(maf_df, reformat_fusion)
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
      "Translation_Start_Site",
      "Nonsense_Mutation",
      "Nonstop_Mutation",
      "In_Frame_Del",
      "In_Frame_Ins",
      "Missense_Mutation",
      "Stop_Codon_Ins",
      "Start_Codon_Del",
      "Fusion",
      "Multi_Hit",
      "Hem_Deletion",
      "Hom_Deletion",
      "Amplification",
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
  removeNonMutated = TRUE,
  annotationFontSize = 0.7,
  SampleNamefontSize = 0.5,
  fontSize = 0.7,
  colors = color_palette
)
dev.off()
