# SNV-caller set up
# 2019
# C. Savonen for ALSF - CCDL
#
# Purpose: Set up for reference files that are used by 01-calculate_vaf_tmb.R
#          script.
#
# Option descriptions
# --annot_rds : File path to where you would like the annotation_rds file to be
#               stored
# --cosmic_og : Path to original COSMIC file. Can be .gz compressed. Give path relative
#               to top directory, 'OpenPBTA-analysis'. Will need to download this from
#               COSMIC at https://cancer.sanger.ac.uk/cosmic/download.
#               These data are available if you register.
# --cosmic_clean : File path specifying where you would like the cleaned brain-related
#                  COSMIC mutations file to be stored.
#
#################################### Set Up ####################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load the optparse library
library(optparse)

# Will need optparse for collecting options
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

################################ Set up options ################################
# Set up optparse options
option_list <- list(
  make_option(
    opt_str = "--cosmic_og", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of 
              'OpenPBTA-analysis') to where the COSMIC mutation file is stored. 
              Will need to download this from COSMIC at 
              https://cancer.sanger.ac.uk/cosmic/download
              These data are available if you register.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--cosmic_clean", type = "character", default = "none",
    help = "Relative file path (assuming from top directory of 
              'OpenPBTA-analysis') where you would like the cleaned COSMIC file
              to be stored.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--annot_rds", type = "character",
    default = "none", help = "Relative file path (assuming from top 
              directory of 'OpenPBTA-analysis') that specifies where you would
              like the annotation_rds file to be stored.",
    metavar = "character"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

##################### Set base directories common file paths ###################
# Make all base names relative to root_dir
opt$cosmic_og <- file.path(root_dir, opt$cosmic_og)
opt$cosmic_clean <- file.path(root_dir, opt$cosmic_clean)
opt$annot_rds <- file.path(root_dir, opt$annot_rds)

################################ Build Annotation ##############################
# We will only run this if it hasn't been run before
if (opt$annot_rds != "none" && !file.exists(opt$annot_rds)) {
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
    opt$annot_rds,
    compress = "gz"
  )
}
######################## Set up COSMIC Mutations file ###########################
# This set up will not be run unless you obtain the original data file
# from https://cancer.sanger.ac.uk/cosmic/download , these data are available
# if you register. The full, unfiltered somatic mutations file
# CosmicMutantExport.tsv for grch38 is used here but only mutations from brain
# related samples are kept.
if (opt$cosmic_clean != "none" && !file.exists(opt$cosmic_clean)) {
  # Print progress message
  message("Setting up COSMIC mutation file. Only need to do this once.")

  # Read in original file
  cosmic_variants <- data.table::fread(opt$cosmic_og,
    data.table = FALSE
  ) %>%
    # Keep only brain mutations so the file is smaller
    dplyr::filter(`Site subtype 1` == "brain") %>%
    # Get rid of spaces in column names
    dplyr::rename_all(dplyr::funs(stringr::str_replace_all(., " ", "_"))) %>%
    # Separate the genome coordinates into their own BED like variables
    dplyr::mutate(
      Chromosome = paste0(
        "chr",
        stringr::word(Mutation_genome_position, sep = ":|-", 1)
      ),
      Start_Position = stringr::word(Mutation_genome_position, sep = ":|-", 2),
      End_Position = stringr::word(Mutation_genome_position, sep = ":|-", 3),
      # Make a base_change variable so we can compare to our set up for PBTA data
      base_change = substr(Mutation_CDS, nchar(Mutation_CDS) - 2, nchar(Mutation_CDS))
    ) %>%
    # Carry over the strand info, but rename to match our PBTA set up
    dplyr::rename(Strand = Mutation_strand) %>%
    # Narrow down to just the needed columns
    dplyr::select(
      "Chromosome", "Start_Position", "End_Position", "Strand",
      "base_change"
    ) %>%
    # Write to a TSV file
    readr::write_tsv(opt$cosmic_clean)
} else {
  warning("A cleaned COSMIC Mutation file was already found with this name. Delete this if you 
         wanted to re-run this step.")
}
