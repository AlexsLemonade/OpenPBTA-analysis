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

# Install GenomicRanges
if (!("GenomicRanges" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}

# Install IRanges
if (!("IRanges" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("Iranges")
}

# Install annotatr 
if (!("annotatr" %in% installed.packages())) {
  install.packages("BiocManager")
  BiocManager::install("annotatr")
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

results_dir <- "results"
if(!dir.exists(results_dir)) {
  dir.create(results_dir)
}

seg_file <-
  data.table::fread(
    file.path(opt$seg_file),
    data.table = FALSE,
    stringsAsFactors = FALSE
  )

#### Format seg file and overlap with hg38 genome annotations ------------------

# Exclude the X and Y chromosomes
seg_no_xy <- subset(seg_file, chrom != "Y" & chrom != "X")
seg_no_xy <- seg_no_xy[, c(2:ncol(seg_no_xy), 1)]
# names(seg_no_xy) <- c("chr", "start", "end", "markers", "CN", "Model")

# Make seg data.frame a GRanges object 
seg_gr <- seg_no_xy %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                          starts.in.df.are.0based = FALSE)

# Define the annotations for the hg38 genome
annotations <- annotatr::build_annotations(genome = "hg38",
                                           annotations = "hg38_basicgenes")

# Create a data.frame with the overlaps between the seg file and hg38 genome 
# annotations 
overlaps <- IRanges::mergeByOverlaps(seg_gr, annotations)

# Remove redundant chr (takes some time, may not be completely necessary)
# overlaps$seg_gr <- gsub("^chr", "", overlaps$seg_gr)
# overlaps$annotations <- gsub("^chr", "", overlaps$annotations)
# overlaps$markers <- NULL

overlaps_df <- as.data.frame(overlaps)

cn_short <- overlaps_df %>%
  dplyr::select(symbol, ID, copy.num) %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::rename(
    "gene_symbol" = "symbol",
    "biospecimen_id" = "ID",
    "copy_number" = "copy.num"
  )

annotated_cn <- cn_short %>%
  dplyr::mutate(label = ifelse(
    copy_number >= 2 & copy_number <= 6,
    "Amplification",
    ifelse(
      copy_number == 0,
      "Hom_Deletion",
      ifelse(copy_number == 1, "Hem_Deletion", copy_number)
    )
  )) %>%
  dplyr::distinct() %>%
  dplyr::select(gene_symbol, biospecimen_id, label, copy_number)

readr::write_tsv(annotated_cn, file.path(results_dir, "annotated_cn.tsv"))
