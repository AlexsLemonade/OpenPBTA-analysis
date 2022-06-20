# This script converts a seg file into a tsv file with CN information and gene
# annotation.
#
# Code adapted from the PPTC PDX Focal Copy NUmber and SV Plots repository here:
# https://github.com/marislab/pptc-pdx-copy-number-and-SVs/blob/master/R/focal-CN-revision.R
#
# Chante Bethell for CCDL 2019 and Jo Lynne Rokita
#
# #### Example Usage
#
# This script is intended to be run via the command line.
# This example assumes it is being run from the root of the repository.
#
# Rscript --vanilla analyses/focal-cn-file-preparation/03-prepare-cn-file.R \
#   --cnv_file data/pbta-cnv-controlfreec.tsv.gz \
#   --gtf_file data/gencode.v27.primary_assembly.annotation.gtf.gz \
#   --metadata data/pbta-histologies.tsv \
#   --filename_lead "controlfreec_annotated_cn"
#   --controlfreec

#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Function ------------------------------------------------------------------

process_annotate_overlaps <- function(cnv_df,
                                      txdb_exons,
                                      filt_na_symbol = TRUE) {
  # This function takes a standardized data.frame that contains genomic range
  # information and finds the overlaps with a TxDb object. It will discard
  # Ensembl gene IDs with no gene symbol by default.
  #
  # Args:
  #   cnv_df: standardized data.frame that contains the segments used in the
  #           CNV caller
  #   txdb_exons: exons to be merged with the CNV data.frame; output of
  #               GenomicFeatures::exons
  #   filt_na_symbol: logical, if TRUE, rows without gene symbols will be
  #                   removed; default is TRUE
  #
  # Returns:
  #   A data.frame with the following columns:
  #     biospecimen_id, status, copy_number, ploidy, ensembl, gene_symbol,
  #     cytoband

  # Make CNV data.frame a GRanges object
  cnv_gr <- cnv_df %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                            starts.in.df.are.0based = FALSE)

  # Create a data.frame with the overlaps between the CNV file and hg38 genome
  # annotations
  overlaps <- IRanges::mergeByOverlaps(cnv_gr, tx_exons)

  # Create a data.frame with selected columns from the `autosome_overlaps`
  # object
  annotated_cn <- data.frame(
    biospecimen_id = overlaps$Kids_First_Biospecimen_ID,
    status = overlaps$status,
    copy_number = overlaps$copy_number,
    ploidy = overlaps$tumor_ploidy,
    ensembl = unlist(overlaps$gene_id),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::distinct() %>%
    # Discard the gene version information in order to get gene symbols and
    # cytoband mappings
    dplyr::mutate(ensembl = gsub("\\..*", "", ensembl))

  annotated_cn <- annotated_cn %>%
    dplyr::mutate(
      gene_symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                          keys = ensembl,
                                          column = "SYMBOL",
                                          keytype = "ENSEMBL"),
      cytoband = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                       keys = ensembl,
                                       column = "MAP",
                                       keytype = "ENSEMBL")
    )

  # drop rows that have NA values for the gene symbol
  if (filt_na_symbol) {
    annotated_cn <- annotated_cn %>%
      dplyr::filter(!is.na(gene_symbol))
  }

  return(annotated_cn)
}

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("--cnv_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains CNV information"
  ),
  optparse::make_option(
    c("--gtf_file"),
    type = "character",
    default = NULL,
    help = "file path to human genome GTF file"
  ),
  optparse::make_option(
    c("--metadata"),
    type = "character",
    default = NULL,
    help = "file path to pbta-histologies.tsv"
  ),
  optparse::make_option(
    c("--filename_lead"),
    type = "character",
    default = "annotated_cn",
    help = "used in file names"
  ),
  optparse::make_option(
    c("--controlfreec"),
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "flag used to indicate if the CNV file is the output of ControlFreeC"
  ),
  optparse::make_option(
    c("--seg"),
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "flag used to indicate if the CNV file was original a SEG file that
            has been prepped by a notebook to include ploidy information"
  ),
  optparse::make_option(
    c("--xy"),
    type = "integer",
    default = 1,
    help = "integer used to indicate if sex chromosome steps should be skipped (0)
            or run (1)"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# error handling related to specifying the CNV method
if (all(opt$controlfreec, opt$seg)) {
  stop("--controlfreec and --seg are mutually exclusive")
}

if (!any(opt$controlfreec, opt$seg)) {
  stop("You must specify the CNV file format by using --controlfreec or
       --seg")
}

# convert xy flag to logical
xy_flag <- as.logical(opt$xy)

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

#### Format CNV file and overlap with hg38 genome annotations ------------------

# we want to standardize the formats between the two methods here and drop
# columns we won't need.
if (opt$seg) {
  cnv_df <- readr::read_tsv(opt$cnv_file) %>%
    dplyr::rename(chr = chrom, start = loc.start, end = loc.end,
                  copy_number = copy.num) %>%
    dplyr::select(-num.mark, -seg.mean) %>%
    dplyr::select(-Kids_First_Biospecimen_ID, dplyr::everything())
}

if (opt$controlfreec) {
  # TODO: filter based on the p-values rather than just dropping them?
  cnv_df <- readr::read_tsv(opt$cnv_file) %>%
    dplyr::rename(copy_number = copy.number) %>%
    dplyr::mutate(chr = paste0("chr", chr)) %>%
    dplyr::select(-segment_genotype, -uncertainty, -WilcoxonRankSumTestPvalue,
                  -KolmogorovSmirnovPvalue) %>%
    dplyr::select(-Kids_First_Biospecimen_ID, dplyr::everything())
}

#### Read in metadata file -----------------------------------------------------

histologies_df <- readr::read_tsv(opt$metadata, guess_max = 10000)

#### Annotation file -----------------------------------------------------------

# this is the output of GenomicFeatures::makeTxDbFromGFF
# TODO: possibly update this when the GTF file gets included in the data
#       download; may also remove the --gtf_file option and hardcode it?
annotation_directory <- file.path(root_dir,
                                  "analyses",
                                  "focal-cn-file-preparation",
                                  "annotation_files")
annotation_file <- file.path(annotation_directory,
                             "txdb_from_gencode.v27.gtf.db")

if (!file.exists(annotation_file)) {
  # Define the annotations for the hg38 genome
  txdb <- GenomicFeatures::makeTxDbFromGFF(
    file = opt$gtf_file,
    format = "gtf"
  )
  # can do this even if the directory exists
  dir.create(annotation_directory, showWarnings = FALSE)
  # write this to file to save time next time
  AnnotationDbi::saveDb(txdb, annotation_file)
} else {
  txdb <- AnnotationDbi::loadDb(annotation_file)
}

# extract the exons but include ensembl gene identifiers
tx_exons <- GenomicFeatures::exons(txdb, columns = "gene_id")

#### Addressing autosomes first ------------------------------------------------

# Exclude the X and Y chromosomes
# Removing copy neutral segments saves on the RAM required to run this step
# and file size
cnv_no_xy <- cnv_df %>%
  dplyr::filter(!(chr %in% c("chrX", "chrY")), status != "neutral")

# Merge and annotated no X&Y
autosome_annotated_cn <- process_annotate_overlaps(cnv_df = cnv_no_xy,
                                                   txdb_exons = tx_exons) %>%
  # mark possible amplifications and deep loss in autosomes
  dplyr::mutate(status = dplyr::case_when(
    copy_number > (2 * ploidy) ~ "amplification",
    copy_number == 0 ~ "deep deletion",
    TRUE ~ status
  ))

# Output file name
autosome_output_file <- paste0(opt$filename_lead, "_autosomes.tsv.gz")

# Save final data.frame to a tsv file
readr::write_tsv(autosome_annotated_cn,
                 file.path(results_dir, autosome_output_file))

#### X&Y -----------------------------------------------------------------------

if (xy_flag) {
  # Filter to just the X and Y chromosomes and remove neutral segments
  # Removing copy neutral segments saves on the RAM required to run this step
  # and file size
  cnv_sex_chrom <- cnv_df %>%
    dplyr::filter(chr %in% c("chrX", "chrY"), status != "neutral")

  # Merge and annotated no X&Y
  sex_chrom_annotated_cn <- process_annotate_overlaps(cnv_df = cnv_sex_chrom,
                                                      txdb_exons = tx_exons) %>%
  # mark possible deep loss in sex chromosome
  dplyr::mutate(status = dplyr::case_when(
    copy_number == 0  ~ "deep deletion",
    TRUE ~ status
  )) 
  

  # Add germline sex estimate into this data.frame
  sex_chrom_annotated_cn <- sex_chrom_annotated_cn %>%
    dplyr::inner_join(dplyr::select(histologies_df,
                                    Kids_First_Biospecimen_ID,
                                    germline_sex_estimate),
                      by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
    dplyr::select(-germline_sex_estimate, dplyr::everything())

  # Output file name
  sex_chrom_output_file <- paste0(opt$filename_lead, "_x_and_y.tsv.gz")

  # Save final data.frame to a tsv file
  readr::write_tsv(sex_chrom_annotated_cn,
                   file.path(results_dir, sex_chrom_output_file))
}
