# Calculate TMB for a given db_file using Strelka and Mutect data
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Option descriptions
#
# --db_file : Path to sqlite database file made from 01-setup_db.py
# --metadata : Relative file path to MAF file to be analyzed. Can be .gz compressed.
#              Assumes file path is given from top directory of 'OpenPBTA-analysis'.
# --coding_regions : File path that specifies the BED regions file that specifies
#                     coding regions that should be used for coding only TMB calculations.
# --overwrite : If specified, will overwrite any files of the same name. Default is FALSE.
# --nonsynfilter_maf: If TRUE, filter out synonymous mutations, keep non-synonymous mutations, based on maftools definition.
# --nonsynfilter_focr: If TRUE, filter out synonymous mutations, keep non-synonymous mutations, based on Friends of Cancer Research Definition.
# --tcga: If TRUE, will skip PBTA metadata specific steps
#
# Command line example:
#
# Rscript analyses/snv-callers/scripts/03-calculate_tmb.R \
# --db_file scratch/testing_snv_db.sqlite \
# --output analyses/snv-callers/results/consensus \
# --metadata data/pbta-histologies.tsv \
# --coding_regions scratch/gencode.v27.primary_assembly.annotation.bed \
# --nonsynfilter_maf
# --overwrite

################################ Initial Set Up ################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "tmb_functions.R"))
source(file.path(root_dir, "analyses", "snv-callers", "util", "split_mnv.R"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load library:
library(optparse)

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-d", "--db_file"), type = "character",
    default = NULL, help = "Path to sqlite database file made from 01-setup_db.py",
    metavar = "character"
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
    opt_str = "--coding_regions", type = "character", default = "none",
    help = "File path that specifies the BED regions file that specifies what
    coding regions should be used for coding only TMB.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--overwrite", action = "store_true",
    default = FALSE, help = "If TRUE, will overwrite any files of
              the same name. Default is FALSE",
    metavar = "character"
  ),
  make_option(
    opt_str = "--nonsynfilter_maf", action = "store_true",
    default = FALSE, help = "If TRUE, filter out synonymous mutations, keep
    non-synonymous mutations, according to maftools definition.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--nonsynfilter_focr", action = "store_true",
    default = FALSE, help = "If TRUE, filter out synonymous mutations, keep
    non-synonymous mutations, according to Friends of Cancer Research definition.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--tcga", action = "store_true",
    default = FALSE, help = "If TRUE, will skip PBTA metadata specific steps",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Make everything relative to root path
opt$metadata <- file.path(root_dir, opt$metadata)
opt$db_file <- file.path(root_dir, opt$db_file)
opt$coding_regions <- file.path(root_dir, opt$coding_regions)

########### Check that the files we need are in the paths specified ############
needed_files <- c(
  opt$metadata, opt$db_file, opt$coding_regions
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

# Get data name
data_name <- ifelse(opt$tcga, "tcga", "pbta")

# Declare output file based on data_name
tmb_coding_file <- file.path(
  opt$output,
  paste0(data_name, "-snv-mutation-tmb-coding.tsv")
)
tmb_all_file <- file.path(
  opt$output,
  paste0(data_name, "-snv-mutation-tmb-all.tsv")
)

# Don't bother if both files exist already and overwrite is FALSE
if (all(file.exists(c(tmb_coding_file, tmb_all_file)), !opt$overwrite)) {
  stop(paste0(tmb_coding_file, tmb_all_file, "both exist and --overwrite is not
              being used. Use --overwrite if you would like to overwrite these files."))
}

######################## Obtain Mutect Strelka mutations #######################
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

# Variant Classification with High/Moderate variant consequences from maftools
maf_nonsynonymous <- c(
  "Missense_Mutation",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Nonsense_Mutation",
  "In_Frame_Del",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)

focr_nonsynonymous <- c(
  "Missense_Mutation",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "Frame_Shift_Ins",
  "Nonsense_Mutation",
  "In_Frame_Del"
)

# Create the consensus for non-MNVs
strelka_mutect_maf_df <- strelka %>%
  # We'll keep the Strelka2 columns and drop Mutect2 columns
  dplyr::inner_join(mutect %>%
    dplyr::select(join_cols),
  by = join_cols,
  copy = TRUE
  ) %>%
  as.data.frame()

# Get Multi-nucleotide calls from mutect as SNVs
split_mutect_df <- split_mnv(mutect) %>%
  dplyr::select(join_cols)

# join MNV calls with strelka
strelka_mutect_mnv <- strelka %>%
  dplyr::inner_join(split_mutect_df,
    by = join_cols,
    copy = TRUE
  ) %>%
  as.data.frame()

if (opt$tcga) {
  strelka_mutect_mnv <- strelka_mutect_mnv %>%
    # In the TCGA MAF files, the Tumor_Sample_Barcode has the biospecimen
    # information but only the first 12 characters are needed to match the metadata
    dplyr::mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 0, 12))
}

# Merge in the MNVs
strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
  dplyr::union(strelka_mutect_mnv,
    by = join_cols
  )

# If the maftools non-synonymous filter is on, filter out synonymous mutations
if (opt$nonsynfilter_maf) {
  strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
  dplyr::filter(Variant_Classification %in% maf_nonsynonymous)
}

# If the FoCR non-synonymous filter is on, filter out synonymous mutations according to that definition
if (opt$nonsynfilter_focr) {
  strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
  dplyr::filter(Variant_Classification %in% focr_nonsynonymous)
}

########################### Set up metadata columns ############################
# Print progress message
message("Setting up metadata...")

# Have to handle TCGA and PBTA metadata differently
if (opt$tcga) {
  # Format two fields of metadata for use with functions
  metadata <- readr::read_tsv(opt$metadata, guess_max = 10000) %>%
    dplyr::mutate(
      short_histology = Primary_diagnosis,
      target_bed_path = file.path(root_dir, "data", BED_In_Use),
      experimental_strategy = "WXS"
    ) %>%
    dplyr::rename(
      Tumor_Sample_Barcode = tumorID,
      target_bed = BED_In_Use
    ) # This field is named differently

  # Manifest files only have first 12 letters of the barcode so we gotta chop the end off
  strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
    dplyr::mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 0, 12))
} else { # pbta data
  # Isolate metadata to only the samples that are in the datasets
  metadata <- readr::read_tsv(opt$metadata, guess_max = 10000) %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% strelka_mutect_maf_df$Tumor_Sample_Barcode) %>%
    dplyr::distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE) %>%
    dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID) %>%
    # Make a Target BED regions column
    dplyr::mutate(
      target_bed = dplyr::recode(experimental_strategy,
        "WGS" = "scratch/intersect_strelka_mutect_WGS.bed",
        "WXS" = "data/WXS.hg38.100bp_padded.bed"
        #TODO: make a padded/unpadded script option
      ),
      target_bed_path = file.path(root_dir, target_bed)
    )

  # Make sure that we have metadata for all these samples.
  if (!all(unique(strelka_mutect_maf_df$Tumor_Sample_Barcode) %in% metadata$Tumor_Sample_Barcode)) {
    stop("There are samples in this MAF file that are not in the metadata.")
  }
}
# Add in metadata
strelka_mutect_maf_df <- strelka_mutect_maf_df %>%
  dplyr::inner_join(metadata %>%
    dplyr::select(
      Tumor_Sample_Barcode,
      experimental_strategy,
      short_histology,
      target_bed,
      target_bed_path
    ),
  by = "Tumor_Sample_Barcode"
  ) %>%
  # Remove samples if they are not WGS or WXS
  dplyr::filter(experimental_strategy %in% c("WGS", "WXS"))

############################# Set Up BED Files #################################
# Make a data.frame of the unique BED file paths and their names
bed_files_key_df <- strelka_mutect_maf_df %>%
  dplyr::select(Tumor_Sample_Barcode, target_bed, target_bed_path) %>%
  dplyr::distinct()

# Get the file paths for the bed files
bed_file_paths <- bed_files_key_df %>%
  dplyr::distinct(target_bed, target_bed_path) %>%
  tibble::deframe()

# Read in each unique BED file and turn into GenomicRanges object
bed_ranges_list <- lapply(bed_file_paths, function(bed_file) {

  # Read in BED file as data.frame
  bed_df <- readr::read_tsv(bed_file,
    col_names = c("chr", "start", "end")
  )

  # Make into a GenomicRanges object
  bed_ranges <- GenomicRanges::GRanges(
    seqnames = bed_df$chr,
    ranges = IRanges::IRanges(
      start = bed_df$start,
      end = bed_df$end
    )
  )
  return(bed_ranges)
})

#################### Set up Coding Region version of BED ranges ################
# Read in the coding regions BED file
coding_regions_df <- readr::read_tsv(opt$coding_regions,
  col_names = c("chr", "start", "end")
)
# Make into a GenomicRanges object
coding_ranges <- GenomicRanges::GRanges(
  seqnames = coding_regions_df$chr,
  ranges = IRanges::IRanges(
    start = coding_regions_df$start,
    end = coding_regions_df$end
  )
)

# For each BED range, find the coding regions intersection
coding_bed_ranges_list <- lapply(bed_ranges_list, function(bed_range,
                                                           coding_grange = coding_ranges) {
  # Find the intersection
  coding_intersect_ranges <- GenomicRanges::intersect(bed_range, coding_grange)

  # Return the reduce version of these ranges
  return(GenomicRanges::reduce(coding_intersect_ranges))
})

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

  # Run TMB calculation on each tumor sample and its respective BED range
  tmb_all_df <- purrr::map2_df(
    bed_files_key_df$Tumor_Sample_Barcode,
    bed_files_key_df$target_bed,
    ~ calculate_tmb(
      tumor_sample_barcode = .x,
      maf_df = strelka_mutect_maf_df,
      bed_ranges = bed_ranges_list[[.y]]
    )
  )

  # Write to TSV file
  readr::write_tsv(tmb_all_df, tmb_all_file)

  # Print out completion message
  message(paste("TMB 'all' calculations saved to:", tmb_all_file))
}
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

  # Run coding TMB calculation on each tumor sample and its
  # respective coding BED range
  tmb_coding_df <- purrr::map2_df(
    bed_files_key_df$Tumor_Sample_Barcode,
    bed_files_key_df$target_bed,
    ~ calculate_tmb(
      tumor_sample_barcode = .x,
      maf_df = strelka_mutect_maf_df,
      bed_ranges = coding_bed_ranges_list[[.y]]
    )
  )

  # Write to TSV file
  readr::write_tsv(tmb_coding_df, tmb_coding_file)

  # Print out completion message
  message(paste("TMB 'coding only' calculations saved to:", tmb_coding_file))
}
