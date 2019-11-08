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

if (!("AnnotationDbi" %in% installed.packages())) {
  BiocManager::install("AnnotationDbi", update = FALSE)
}

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

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
    c("--metadata"),
    type = "character",
    default = NULL,
    help = "file path to file that contains the sample metadata"
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
    c("--cnvkit"),
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "flag used to indicate if the CNV file is the output of CNVkit"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# error handling related to specifying the CNV method
if (all(opt$controlfreec, opt$cnvkit)) {
  stop("--controlfreec and --cnvkit are mutually exclusive")
}

if (!any(opt$controlfreec, opt$cnvkit)) {
  stop("You must specify the CNV file format by using --controlfreec or
       --cnvkit")
}

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
# columns we won't need to
if (opt$cnvkit) {
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

# Remove neutral copy number calls and mark possible amplifcation
cnv_df <- cnv_df %>%
  dplyr::filter(status != "neutral") %>%
  dplyr::mutate(status = dplyr::case_when(
    copy_number > (2 * tumor_ploidy) ~ "amplification",
    TRUE ~ status
  ))

# Addressing autosomes first
# Exclude the X and Y chromosomes
cnv_no_xy <- cnv_df %>%
  dplyr::filter(!(chr %in% c("chrX", "chrY")))

# Make seg data.frame a GRanges object
cnv_no_xy_gr <- cnv_no_xy %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                                          starts.in.df.are.0based = FALSE)

# Define the annotations for the hg38 genome
txdb <- GenomicFeatures::makeTxDBFromGFF(
  file = opt$gtf_file,
  format = "gtf"
)

# we'll only look at genes on chromosomes 1:22 -- this will get rid of things
# like alternates
chroms <- paste0("chr", 1:22)
chrom_filter <- list(tx_chrom = chroms)

# extract the genes using the chromosome filter
chr_exons <- GenomicFeatures::exons(txdb,  filter = chrom_filter)

# Create a data.frame with the overlaps between the seg file and hg38 genome
# annotations
overlaps <- IRanges::mergeByOverlaps(cnv_no_xy_gr, chr_exons)

overlap_symbols <-
  AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = overlaps$gene_id,
    column = "SYMBOL",
    keytype = "ENTREZID"
  )

overlaps_cytoband <-
  AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = overlaps$gene_id,
    column = "MAP",
    keytype = "ENTREZID"
  )

overlaps$gene_id <- overlap_symbols
overlaps$cytoband <- overlaps_cytoband

# Create a data.frame with selected columns from the `overlaps` object
cn_short <- data.frame(
  gene_symbol = overlaps$gene_id,
  biospecimen_id = overlaps$ID,
  copy_number = overlaps$copy.num) %>%
  dplyr::filter(!(is.na(gene_symbol)))

# Create a `label` column with values specifying the type of CN aberration and
# save as a data.frame with `gene_symbol`, `biospecimen_id`, and `copy_number`
annotated_cn <- cn_short %>%
  dplyr::distinct() %>%
  dplyr::mutate(label = dplyr::case_when(
    copy_number == 3 | copy_number == 4  ~ "Gain",
    copy_number == 0 ~ "Hom_Deletion",
    copy_number == 1 ~ "Hem_Deletion",
    copy_number >= 5 ~ "Amplification")) %>%
  dplyr::select(gene_symbol, biospecimen_id, label, copy_number)

# Save final data.frame to a tsv file
readr::write_tsv(annotated_cn, file.path(results_dir, "annotated_cn.tsv"))
