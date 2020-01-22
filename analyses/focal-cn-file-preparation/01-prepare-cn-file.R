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

#### Function ------------------------------------------------------------------

process_annotate_overlaps <- function(cnv_df,
                                      cds_gr,
                                      filt_na_symbol = TRUE) {
  # This function takes a standardized data.frame that contains genomic range
  # information and finds the overlaps with a CDS GRanges object. It will
  # discard Ensembl gene IDs with no gene symbol by default.
  #
  # Args:
  #   cnv_df: standardized data.frame that contains the segments used in the
  #           CNV caller
  #   cds_gr: `GenomicRanges` object containing coding sequences to be merged
  #           with the CNV data.frame; output of `run-prepare-cn.sh`
  #   filt_na_symbol: logical, if TRUE, rows without gene symbols will be
  #                   removed; default is TRUE
  #
  # Returns:
  #   A data.frame with the following columns:
  #     biospecimen_id, status, copy_number, ploidy, ensembl, gene_symbol,
  #     cytoband

  # Make CNV data.frame a GRanges object
  cnv_gr <- cnv_df %>%
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      starts.in.df.are.0based = FALSE
    )

  # Find the overlaps between the CNV file and GENCODE v27 filtered hg38
  # annotations file (converted to bed file in `run-prepare-cn.sh`).
  overlaps <- GenomicRanges::findOverlaps(cnv_gr, cds_gr)

  # Carry over gene ids to `cnv_gr` object
  cnv_gr$gene_id <- rep(NA, length(cnv_gr))

  cnv_gr$gene_id[overlaps@from] <- cds_gr$gene_id[overlaps@to]

  # Find the overlaps between the CNV file and UCSC cytoband data
  overlaps_ucsc <- GenomicRanges::findOverlaps(cnv_gr, ucsc_cytoband_gr)

  # Create a data.frame with selected columns from the `overlaps`
  # object
  annotated_cn <- data.frame(
    biospecimen_id = cnv_gr$Kids_First_Biospecimen_ID[overlaps_ucsc@from],
    status = cnv_gr$status[overlaps_ucsc@from],
    copy_number = cnv_gr$copy_number[overlaps_ucsc@from],
    ploidy = cnv_gr$tumor_ploidy[overlaps_ucsc@from],
    ensembl = cnv_gr$gene_id[overlaps_ucsc@from],
    cytoband = ucsc_cytoband_gr@elementMetadata@listData$cytoband[overlaps_ucsc@to],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!(is.na(ensembl))) %>%
    # Discard the gene version information in order to get gene symbols and
    # cytoband mappings
    dplyr::mutate(ensembl = gsub("\\..*", "", ensembl))

  annotated_cn <- annotated_cn %>%
    dplyr::mutate(
      gene_symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
        keys = ensembl,
        column = "SYMBOL",
        keytype = "ENSEMBL"
      )
    ) %>%
    dplyr::select(-cytoband) %>%
    dplyr::distinct()

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
    c("--cds_file"),
    type = "character",
    default = NULL,
    help = "file path to human genome cds filtered GTF file"
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
    c("--cnvkit"),
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "flag used to indicate if the CNV file is the output of CNVkit"
  ),
  optparse::make_option(
    c("--gistic"),
    type = "logical",
    action = "store_true",
    default = FALSE,
    help = "flag used to indicate if the GISTIC file should be incorporated"
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
if (all(opt$controlfreec, opt$cnvkit)) {
  stop("--controlfreec and --cnvkit are mutually exclusive")
}

if (!any(opt$controlfreec, opt$cnvkit)) {
  stop("You must specify the CNV file format by using --controlfreec or
       --cnvkit")
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

if (opt$gistic) {
  # Path to gistic zip file results
  gistic_zip <-
    file.path(root_dir, "data", "pbta-cnv-cnvkit-gistic.zip")
  
  # Get the name of the folder first since this changes with updates to the date
  folder_name <- unzip(zipfile = gistic_zip, list = TRUE) %>%
    dplyr::filter(Length == 0) %>%
    dplyr::pull("Name")
  
  # Path to gistic results directory
  gistic_results <-
    file.path(root_dir,
              "analyses",
              "focal-cn-file-preparation",
              "gistic-results")
  
  # Extract files if it hasn't been done yet
  if (!dir.exists(gistic_results)) {
    # Unzip but only extract the files not the folder
    unzip(zipfile = gistic_zip,
          exdir = gistic_results)
  }
  # Designate GISTIC results folder to extract to
  gistic_results <- file.path(gistic_results, folder_name)
  
  # Read in broad values by arm GISTIC output file
  gistic_df <-
    data.table::fread(
      paste0(
        gistic_results,
        "broad_values_by_arm.txt"
      ),
      data.table = FALSE
    )
}

# Read in UCSC cytoband data. The decision to implement the UCSC hg38 cytoband
# file was made based on a comparison done between the cytoband calls in the
# `org.Hs.eg.db` package and the calls in the UCSC file. We found that they
# disagreed in ~11,800 calls out of ~800,000 and the `UCSC file` contains more
# cytoband calls.
ucsc_cytoband <-
  data.table::fread(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz"
  )

# Make the UCSC cytoband data.frame a GRanges object
ucsc_cytoband_gr <- ucsc_cytoband %>%
  dplyr::select(chr = V1, start = V2, end = V3, cytoband = V4) %>%
  dplyr::mutate(cytoband = paste0(gsub("chr", "", chr), cytoband)) %>%
  GenomicRanges::makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = FALSE
  )

#### Format CNV file and overlap with hg38 genome annotations ------------------

# We want to standardize the formats between the two methods here and drop
# columns we won't need to
if (opt$cnvkit) {
  cnv_df <- readr::read_tsv(opt$cnv_file) %>%
    dplyr::rename(
      chr = chrom, start = loc.start, end = loc.end,
      copy_number = copy.num
    ) %>%
    dplyr::select(-num.mark, -seg.mean) %>%
    dplyr::select(-Kids_First_Biospecimen_ID, dplyr::everything())
}

if (opt$controlfreec) {
  # TODO: filter based on the p-values rather than just dropping them?
  cnv_df <- readr::read_tsv(opt$cnv_file) %>%
    dplyr::rename(copy_number = copy.number) %>%
    dplyr::mutate(chr = paste0("chr", chr)) %>%
    dplyr::select(
      -segment_genotype, -uncertainty, -WilcoxonRankSumTestPvalue,
      -KolmogorovSmirnovPvalue
    ) %>%
    dplyr::select(-Kids_First_Biospecimen_ID, dplyr::everything())
}

#### Read in metadata file -----------------------------------------------------

histologies_df <- readr::read_tsv(opt$metadata)

#### Annotation file -----------------------------------------------------------

# Read in prepared cds file
cds_df <- data.table::fread(opt$cds_file, data.table = FALSE)

# Create the GRanges object
cds_gr <- cds_df %>%
  dplyr::select(
    seqnames = V1,
    start = V2,
    end = V3,
    gene_id = V4
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = FALSE
  )

#### Addressing autosomes first ------------------------------------------------

# Exclude the X and Y chromosomes
cnv_no_xy <- cnv_df %>%
  dplyr::filter(!(chr %in% c("chrX", "chrY")), status != "neutral")

# Merge and annotated no X&Y
autosome_annotated_cn <- process_annotate_overlaps(
  cnv_df = cnv_no_xy,
  cds_gr = cds_gr
) %>%
  # mark possible amplifications in autosomes
  dplyr::mutate(status = dplyr::case_when(
    copy_number > (2 * ploidy) ~ "amplification",
    TRUE ~ status
  ))

#### Wrangle GISTIC data for autosome data.frame-------------------------------

if (opt$gistic) {
  # Transpose the GISTIC data.frame
  transposed_gistic_df <- gistic_df %>%
    tibble::column_to_rownames("Chromosome Arm") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("biospecimen_id") %>%
    # Gather the chromosome arm data into one column and the GISTIC call in another
    tidyr::gather(key = "arm", value = "gistic_call", -biospecimen_id) %>%
    # Filter the GISTIC data to only the overlapping samples
    dplyr::filter(biospecimen_id %in% autosome_annotated_cn$biospecimen_id) %>%
    # Recode the gistic
    dplyr::mutate(broad_status = dplyr::case_when(
      gistic_call == -1 ~ "loss",
      gistic_call == 0 ~ "neutral",
      gistic_call == 1 ~ "gain",
      TRUE ~ "amplification"
    )) %>%
    dplyr::select(-gistic_call)

  # Add GISTIC data to final data.frame and remove sex chromosomes
  autosome_annotated_cn <- transposed_gistic_df %>%
    dplyr::filter(!(arm %in% c("Xp", "Xq", "Yp", "Yq"))) %>%
    dplyr::inner_join(autosome_annotated_cn, by = "biospecimen_id")
}

# Output file name
autosome_output_file <- paste0(opt$filename_lead, "_autosomes.tsv.bz2")

# Save final data.frame to a tsv file
readr::write_tsv(
  autosome_annotated_cn,
  file.path(results_dir, autosome_output_file)
)

#### X&Y -----------------------------------------------------------------------

if (xy_flag) {
  # Filter to just the X and Y chromosomes
  cnv_sex_chrom <- cnv_df %>%
    dplyr::filter(chr %in% c("chrX", "chrY"))

  # Merge and annotated no X&Y
  sex_chrom_annotated_cn <- process_annotate_overlaps(
    cnv_df = cnv_sex_chrom,
    cds_gr = cds_gr
  )

  # Add germline sex estimate into this data.frame
  sex_chrom_annotated_cn <- sex_chrom_annotated_cn %>%
    dplyr::inner_join(dplyr::select(
      histologies_df,
      Kids_First_Biospecimen_ID,
      germline_sex_estimate
    ),
    by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")
    ) %>%
    dplyr::select(-germline_sex_estimate, dplyr::everything())

  if (opt$gistic) {
    #### Wrangle GISTIC data for sex chromosome data.frame-----------------------

    # Add GISTIC data to final data.frame and select only the sex chromsomes
    sex_chrom_annotated_cn <- transposed_gistic_df %>%
      dplyr::filter(arm %in% c("Xp", "Xq", "Yp", "Yq")) %>%
      dplyr::inner_join(sex_chrom_annotated_cn, by = "biospecimen_id")
  }

  # Output file name
  sex_chrom_output_file <- paste0(opt$filename_lead, "_x_and_y.tsv.bz2")

  # Save final data.frame to a tsv file
  readr::write_tsv(
    sex_chrom_annotated_cn,
    file.path(results_dir, sex_chrom_output_file)
  )
}
