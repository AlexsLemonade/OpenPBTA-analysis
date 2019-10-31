# This script converts a seg file into a tsv file with CN information and gene
# annotation.
#
# Code adapted from the PPTC PDX Focal Copy NUmber and SV Plots repository here:
# https://github.com/marislab/pptc-pdx-copy-number-and-SVs/blob/master/R/focal-CN-revision.R
#
# Chante Bethell for CCDL 2019 and Jo Lynne Rokita
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/oncoprint-landscape/01-prepare-cn-file.R'

#### Set Up --------------------------------------------------------------------

# We require bioconductor
if (!("BiocManager" %in% installed.packages())) {
  install.packages("BiocManager")
}

# Install GenomicRanges
if (!("GenomicRanges" %in% installed.packages())) {
  BiocManager::install("GenomicRanges", update = FALSE)
}

# Install IRanges
if (!("IRanges" %in% installed.packages())) {
  BiocManager::install("IRanges", update = FALSE)
}

# Install annotatr 
if (!("annotatr" %in% installed.packages())) {
  BiocManager::install("annotatr", update = FALSE)
}

# hg38 genome annotations
if (!("TxDb.Hsapiens.UCSC.hg38.knownGene" %in% installed.packages())) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", update = FALSE)
}

if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-f", "--seg_file"),
    type = "character",
    default = NULL,
    help = "file path to SEG file that contains cnv information"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results directory 
results_dir <-
  file.path(root_dir, "analyses", "focal-cn-file-preparation", "results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read in data from seg file
seg_df <-
  data.table::fread(
    opt$seg_file,
    data.table = FALSE,
    stringsAsFactors = FALSE
  )

#### Format seg file and overlap with hg38 genome annotations ------------------

# Exclude the X and Y chromosomes and rearrange ID column to be the last column
seg_no_xy <- seg_df %>%
  dplyr::filter(!(chrom %in% c("chrX", "chrY"))) %>%
  dplyr::select(-ID, dplyr::everything())

# Make seg data.frame a GRanges object 
seg_gr <- seg_no_xy %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                          starts.in.df.are.0based = FALSE)

# Define the annotations for the hg38 genome
annotations <- annotatr::build_annotations(genome = "hg38",
                                           annotations = "hg38_genes_exons")

# Create a data.frame with the overlaps between the seg file and hg38 genome 
# annotations 
overlaps <- IRanges::mergeByOverlaps(seg_gr, annotations)

overlaps_df <- as.data.frame(overlaps)

cn_short <- overlaps_df %>%
  dplyr::select(symbol, ID, copy.num) %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::rename(
    "gene_symbol" = "symbol",
    "biospecimen_id" = "ID",
    "copy_number" = "copy.num"
  )

# Create a `label` column with values specifying the type of CN aberration and
# save as a data.frame with `gene_symbol`, `biospecimen_id`, and `copy_number`
annotated_cn <- cn_short %>%
  dplyr::distinct() %>%
  dplyr::mutate(label = dplyr::case_when(
    copy_number >= 2 & copy_number <= 6  ~ "Amplification",
    copy_number == 0 ~ "Hom_Deletion",
    copy_number == 1 ~ "Hem_Deletion")) %>%
  dplyr::select(gene_symbol, biospecimen_id, label, copy_number)

# Save final data.frame to a tsv file
readr::write_tsv(annotated_cn, file.path(results_dir, "annotated_cn.tsv"))
