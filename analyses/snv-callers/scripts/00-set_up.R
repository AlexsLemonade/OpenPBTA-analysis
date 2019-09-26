# SNV Caller set up and functions
#
# 2019
# C. Savonen for ALSF - CCDL
#
# Purpose: Set up for analyzing MAF files
#
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "wrangle_functions.R"))
source(file.path(root_dir, "analyses", "snv-callers", "util", "plot_functions.R"))

#################################### Set Up ####################################

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Will need R.utils for zipping up the results file
if (!("optparse" %in% installed.packages())) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
}
# Will need R.utils for zipping up the results file
if (!("R.utils" %in% installed.packages())) {
  install.packages("R.utils", repos = "http://cran.us.r-project.org")
}
# Install package if not installed
if (!("annotatr" %in% installed.packages())) {
  BiocManager::install("annotatr")
}
# Install package if not installed
if (!("data.table" %in% installed.packages())) {
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}

##################### Set base directories common file paths ###################
# Declare base directory names
scratch_dir <- file.path(root_dir, "scratch")
data_dir <- file.path(root_dir, "data")
snv_dir <- file.path(root_dir, "analyses", "snv-callers")

# Directories that we will need
bed_dir <- file.path(snv_dir, "bed_regions")
base_results_dir <- file.path(snv_dir, "results")
base_plots_dir <- file.path(snv_dir, "plots")
cosmic_dir <- file.path(snv_dir, "cosmic")

# Create these folders if they haven't been created yet
if (!dir.exists(bed_dir)) {
  dir.create(bed_dir)
}
if (!dir.exists(base_results_dir)) {
  dir.create(base_results_dir)
}
if (!dir.exists(base_plots_dir)) {
  dir.create(base_plots_dir)
}
if (!dir.exists(cosmic_dir)) {
  dir.create(cosmic_dir)
}

# Declare reference file names
original_metadata <- file.path(data_dir, "pbta-histologies.tsv")

# These data were obtained from https://cancer.sanger.ac.uk/cosmic/download
# These data are available if you register.
# The full, unfiltered somatic mutations file CosmicMutantExport.tsv.gz for grch38
# is used here.
cosmic_file <- file.path(snv_dir, "CosmicMutantExport.tsv.gz")

# This is the cleaned version that only contains the genomic coordinates and
# base changes from the original file for brain sample mutations only.
cosmic_clean_file <- file.path(cosmic_dir, "brain_cosmic_variants_coordinates.tsv")

# Will set up this annotation file below
annot_rds <- file.path(scratch_dir, "hg38_genomic_region_annotation.rds")

# This will be used for all samples, but WGS bed regions are caller-specific
wxs_bed_file <- file.path(bed_dir, "wxs_bed_regions_all.tsv")

######################## Set up WXS BED regions file ###########################
# WGS BED regions files are caller specific so will be set up separately
if (!file.exists(wxs_bed_file)) {
  # This is a link to the exome BED file provided by
  wxs_bed_url <- "https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/Exome-AZ_V2.bed"

  # Download the file with these regions
  download.file(wxs_bed_url, destfile = wxs_bed_file)
}

################################ Build Annotation ##############################
# We will only run this if it hasn't been run before
if (!file.exists(annot_rds)) {
  # Print progress message
  message("Setting up genomic region annotation file. Only need to do this once.")

  # Here we are specifying what types of annotation we would like
  annots <- c(
    "hg38_basicgenes",
    "hg38_genes_intergenic",
    "hg38_genes_intronexonboundaries",
    "hg38_genes_exonintronboundaries",
    "hg38_genes_5UTRs",
    "hg38_genes_exons",
    "hg38_genes_introns",
    "hg38_genes_3UTRs"
  )

  # Build all the annotations into GRanges object
  annotations <- annotatr::build_annotations(genome = "hg38", annotations = annots)

  # Write this object so we don't have to write it again
  readr::write_rds(annotations,
    annot_rds,
    compress = "gz"
  )
}

######################## Set up COSMIC Mutations file ###########################
# This set up will not be run unless you obtain the original data file
# from https://cancer.sanger.ac.uk/cosmic/download , these data are available
# if you register. The full, unfiltered somatic mutations file
# CosmicMutantExport.tsv for grch38 is used here but only mutations from brain
# related samples are kept.
if (!file.exists(cosmic_clean_file)) {
  # Print progress message
  message("Setting up COSMIC mutation file. Only need to do this once.")

  # Read in original file
  cosmic_variants <- data.table::fread(cosmic_file,
    data.table = FALSE
  ) %>%
    # Keep only brain mutations so the file is smaller
    dplyr::filter(`Site subtype 1` == "brain")
  # Get rid of spaces in column names
  cosmic_variants <- cosmic_variants %>%
    dplyr::rename_all(dplyr::funs(stringr::str_replace_all(., " ", "_"))) %>%
    # Separate the genome coordinates into their own BED like variables
    dplyr::mutate(
      Chromosome = paste0(
        "chr",
        stringr::word(Mutation_genome_position, sep = ":", 1)
      ),
      Start_Position = stringr::word(Mutation_genome_position, sep = ":|-", 1),
      End_Position = stringr::word(Mutation_genome_position, sep = "-", 2),
      # Make a base_change variable so we can compare to our set up for PBTA data
      base_change = substr(Mutation_CDS, nchar(Mutation_CDS) - 2, 10)
    ) %>%
    # Carry over the Strand info,  but rename
    dplyr::rename(Strand = Mutation_strand) %>%
    # Narrow down to just the needed columns
    dplyr::select(
      "Chromosome", "Start_Position", "End_Position", "Strand",
      "base_change"
    ) %>%
    # Write to a TSV file
    readr::write_tsv(cosmic_clean_file)
}

################### Set up caller-specific WGS BED regions files ###############
# **NOTE:** this section cannot be completed until we receive the WGS files for
# the other callers.

# Set up strelka2 WGS bed regions file
if (!file.exists(file.path(bed_dir, "strelka2_wgs_bed_regions.tsv"))) {
  # Get the strelka2 regions
  tmp <- readLines("https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/3#issuecomment-530583149")

  # Write them to a file.
  text <- grep("^chr", tmp, value = TRUE)
  text <- substr(text, 1, nchar(text) - 1)

  data.frame(
    stringr::word(text, sep = "\t"),
    stringr::word(text, sep = "\t", 2),
    stringr::word(text, sep = "\t", 3)
  ) %>%
    readr::write_tsv(file.path(bed_dir, "strelka2_wgs_bed_regions.tsv"),
      col_names = FALSE
    )

  # Get rid of temporary object
  rm(tmp)
}

# Set up mutect2 WGS bed regions file
if (!file.exists(file.path(bed_dir, "mutect2_wgs_bed_regions.tsv"))) {
  # This is a link to the BED regions file from @migbro.
  bed_url <- "https://github.com/AlexsLemonade/OpenPBTA-analysis/files/3602922/wgs_canonical_calling_regions.hg38.txt"

  # Download the file with these regions
  download.file(bed_url, destfile = file.path(bed_dir, "mutect2_wgs_bed_regions.tsv"))
}
