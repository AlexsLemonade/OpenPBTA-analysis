# Calculate TMB for a given MAF file
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# --consensus : File path to the MAF-like file.
# --db_file : Path to sqlite database file made from 01-setup_db.py
# --gtf_file : File path to human genome GTF file. Used to calculate coding genome size.
# --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
#              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --bed_wgs : File path that specifies the caller-specific
#             BED regions file. If two files are given, separate their file paths with a
#             comma. The intersection of these ranges will be taken and used as a denominator. 
#             Assumes from top directory, 'OpenPBTA-analysis'
# --bed_wxs : File path that specifies the WXS BED regions file. Assumes file path
#             is given from top directory of 'OpenPBTA-analysis'
# --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
#
# Command line example:
#
# Rscript analyses/snv-callers/scripts/03-calculate_tmb.R \
#   --consensus analyses/snv-callers/results/consensus/consensus_snv.maf.tsv \
#   --db_file $dbfile \
#   --gtf_file data/gencode.v27.primary_assembly.annotation.gtf.gz \
#   --output analyses/snv-callers/results/consensus \
#   --metadata data/pbta-histologies.tsv \
#   --bed_wgs data/WGS.hg38.strelka2.unpadded.bed,data/WGS.hg38.mutect2.unpadded.bed  \
#   --bed_wxs data/WXS.hg38.100bp_padded.bed \
#   --overwrite

################################ Initial Set Up ################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "wrangle_functions.R"))
source(file.path(root_dir, "analyses", "snv-callers", "util", "split_mnv.R"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-c", "--consensus"), type = "character",
    default = NULL, help = "File path to the MAF-like file",
    metavar = "character"
  ),
  make_option(
    opt_str = c("-d", "--db_file"), type = "character",
    default = NULL, help = "Path to sqlite database file made from 01-setup_db.py",
    metavar = "character"
  ),
  optparse::make_option(
    c("--gtf_file"),
    type = "character",
    default = NULL,
    help = "File path to human genome GTF file. Used to calculate coding genome size."
  ),
  make_option(
    opt_str = c("-o", "--output"), type = "character",
    default = NULL, help = "Path to folder where you would like the
              output TMB file from this script to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--metadata", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of
              'OpenPBTA-analysis') to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--bed_wgs", type = "character", default = "none",
    help = "File path that specifies the caller-specific
    BED regions file. If two files are given, separate their file paths with a
    comma. The intersection of these ranges will be taken and used as a denominator. 
    Assumes from top directory, 'OpenPBTA-analysis'",
    metavar = "character"
  ),
  make_option(
    opt_str = "--bed_wxs", type = "character", default = "none",
    help = "File path that specifies the WXS BED regions file. Assumes
              from top directory, 'OpenPBTA-analysis'",
    metavar = "character"
  ),
  make_option(
    opt_str = "--overwrite", action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any files of
              the same name. Default is FALSE",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Split anywhere there are commas for multiple BED files
opt$bed_wgs <- unlist(strsplit(opt$bed_wgs, ","))

# Make everything relative to root path
opt$consensus <- file.path(root_dir, opt$consensus)
opt$db_file <- file.path(root_dir, opt$db_file)
opt$gtf_file <- file.path(root_dir, opt$gtf_file)
opt$metadata <- file.path(root_dir, opt$metadata)
opt$bed_wgs <- file.path(root_dir, opt$bed_wgs)
opt$bed_wxs <- file.path(root_dir, opt$bed_wxs)

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$consensus, opt$metadata, opt$bed_wgs, opt$bed_wxs, opt$db_file, opt$gtf_file
)

# Get list of which files were found
files_found <- file.exists(needed_files)

# Report error if any of them aren't found
if (!all(files_found)) {
  stop(paste("\n Could not find needed file(s):",
    needed_files[which(!files_found)],
    "Check your options and set up.",
    sep = "\n"
  ))
}

############################### Set Up Output #####################################
# Set and make the plots directory
opt$output <- file.path(root_dir, opt$output)

# Make output folder
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Declare output file paths
tmb_coding_file <- file.path(opt$output, "consensus_snv_tmb_coding_only.tsv")
tmb_all_file <- file.path(opt$output, "consensus_snv_tmb_all.tsv")

# Don't bother if both files exist already and overwrite is FALSE
if (all(file.exists(c(tmb_coding_file, tmb_all_file)), !opt$overwrite)) {
  stop(paste0(tmb_coding_file, tmb_all_file, "both exist and --overwrite is not
              being used. Use --overwrite if you would like to overwrite these files."))
}
########################### Set up this consensus data ##########################
# Print progress message
message(paste("Reading in", opt$consensus, "MAF data..."))

# Read in this MAF, skip the version number
maf_df <- data.table::fread(opt$consensus, data.table = FALSE)

# Print progress message
message("Setting up metadata...")

# Isolate metadata to only the samples that are in the datasets
metadata <- readr::read_tsv(opt$metadata) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% maf_df$Tumor_Sample_Barcode) %>%
  dplyr::distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE) %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID)

# Make sure that we have metadata for all these samples.
if (!all(unique(maf_df$Tumor_Sample_Barcode) %in% metadata$Tumor_Sample_Barcode)) {
  stop("There are samples in this MAF file that are not in the metadata.")
}

# Add the experimental strategy column on the data.frame for calculating purposes
maf_df <- maf_df %>%
  dplyr::inner_join(metadata %>%
    dplyr::select(
      Tumor_Sample_Barcode,
      experimental_strategy,
      short_histology
    ))

############################ Set up this exon regions ##########################
# Declare file path
annotation_file <- file.path("scratch", "txdb_from_gencode.v27.gtf.db")

# Only remake the file if it doesn't exist
if (!file.exists(annotation_file)) {
  message("Creating exon annotation file")

  # Define the annotations for the hg38 genome
  txdb <- GenomicFeatures::makeTxDbFromGFF(
    file = opt$gtf_file,
    format = "gtf"
  )

  # Write this to file to save time next time
  AnnotationDbi::saveDb(txdb, annotation_file)
} else {
  txdb <- AnnotationDbi::loadDb(annotation_file)
}

# extract the exons but include ensembl gene identifiers
tx_exons <- GenomicFeatures::exons(txdb, columns = "gene_id")

# Extract the ranges and sum to get the size
coding_genome_size <- sum(GenomicRanges::width(
  GenomicRanges::reduce(
    tx_exons
  )
))

############################# Calculate TMB ####################################

############################# Coding TMB file ##################################
# If the file exists or the overwrite option is not being used, run TMB calculations
if (file.exists(tmb_coding_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The 'coding only' Tumor Mutation Burden file already exists: \n",
    tmb_coding_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Print out warning if this file is going to be overwritten
  if (file.exists(tmb_coding_file)) {
    warning("Overwriting existing 'coding only' TMB file.")
  }

  # Print out progress message
  message(paste("Calculating 'coding only' TMB..."))

  # Filter out mutations that are outside of these regions.
  coding_maf_df <- snv_ranges_filter(maf_df, keep_ranges = tx_exons)

  # Calculate coding only TMBs and write to file
  tmb_coding_df <- calculate_tmb(coding_maf_df,
    wgs_size = coding_genome_size,
    wxs_size = coding_genome_size
  )
  readr::write_tsv(tmb_coding_df, tmb_coding_file)

  # Print out completion message
  message(paste("TMB 'coding only' calculations saved to:", tmb_coding_file))
}

########################### All mutations TMB file #############################
# If the file exists or the overwrite option is not being used, run TMB calculations
if (file.exists(tmb_all_file) && !opt$overwrite) {
  # Stop if this file exists and overwrite is set to FALSE
  warning(cat(
    "The 'all mutations' Tumor Mutation Burden file already exists: \n",
    tmb_all_file, "\n",
    "Use --overwrite if you want to overwrite it."
  ))
} else {
  # Print out warning if this file is going to be overwritten
  if (file.exists(tmb_coding_file)) {
    warning("Overwriting existing 'all mutations' TMB file.")
  }
  #### Calculate intersection genome size
  # Read in BED region files for TMB calculations
  wgs_beds <- lapply(opt$bed_wgs, function(bed) {

    # Read in BED formatted file
    bed <- readr::read_tsv(bed, col_names = FALSE)

    # Create the GRanges object from bed
    GenomicRanges::GRanges(
      seqnames = bed$X1,
      ranges = IRanges::IRanges(
        start = bed$X2,
        end = bed$X3
      )
    )
  })
  # Get intersection
  intersection_ranges <- GenomicRanges::intersect(wgs_beds[[1]], wgs_beds[[2]])

  # Get genome size
  intersect_genome_size <- sum(GenomicRanges::width(
    GenomicRanges::reduce(
      intersection_ranges
    )
  ))

  ######## Obtain Mutect Strelka intersection
  # Start up connection
  con <- DBI::dbConnect(RSQLite::SQLite(), opt$db_file)

  # Designate caller tables from SQL file
  strelka <- dplyr::tbl(con, "strelka")
  mutect <- dplyr::tbl(con, "mutect")

  # Specify the columns to join by
  join_cols <- c(
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Allele",
    "Tumor_Sample_Barcode"
  )

  # Create the consensus for non-MNVs
  strelka_mutect_maf_df <- strelka %>%
    dplyr::inner_join(mutect, by = join_cols)

  # Get Multi-nucleotide calls from mutect as SNVs
  split_mutect_df <- split_mnv(mutect)

  # join MNV calls with strelka
  strelka_mutect_mnv <- strelka %>%
    dplyr::inner_join(split_mutect_df,
      by = join_cols,
      copy = TRUE
    ) %>%
    as.data.frame()

  # Add in the MNVs
  strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
    dplyr::union_all(strelka_mutect_mnv,
      copy = TRUE
    ) %>%
    dplyr::inner_join(metadata %>%
      dplyr::select(
        "Tumor_Sample_Barcode",
        "experimental_strategy",
        "short_histology"
      ),
    copy = TRUE
    ) %>%
    as.data.frame()

  # .x is messing up the maf_to_granges function
  colnames(strelka_mutect_maf_df) <- gsub("\\.x$", "", colnames(strelka_mutect_maf_df))

  # For WXS samples, filter out mutations that are outside of these coding regions.
  filt_wxs_maf_df <- snv_ranges_filter(dplyr::filter(
    strelka_mutect_maf_df,
    experimental_strategy == "WXS"
  ),
  keep_ranges = tx_exons
  )

  # Bind the filtered WXS sample rows back to the WGS samples
  strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
    dplyr::filter(experimental_strategy == "WGS") %>% # Note that `Panel` samples are not included here. 
    dplyr::bind_rows(filt_wxs_maf_df)

  # Calculate TMBs and write to TMB file
  tmb_all_df <- calculate_tmb(strelka_mutect_maf_df,
    wgs_size = intersect_genome_size,
    wxs_size = coding_genome_size
  )
  readr::write_tsv(tmb_all_df, tmb_all_file)

  # Print out completion message
  message(paste("TMB 'all' calculations saved to:", tmb_all_file))
}
