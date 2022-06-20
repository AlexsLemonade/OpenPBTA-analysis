# This script subsets the focal copy number, RNA expression, histologies`
# files to include only Chordoma samples.

# Written originally Chante Bethell 2019
# (Adapted for this module by Candace Savonen 2020)
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/molecular-subtyping-chordoma/00-subset-files-for-chordoma.R'


#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results and plots directory
results_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-chordoma", "chordoma-subset")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies-base.tsv"), guess_max = 10000)

#### Filter metadata -----------------------------------------------------------
# Select wanted columns in metadata for merging and assign to a new object
chordoma_metadata <- metadata %>%
  dplyr::filter(pathology_diagnosis == "Chordoma") %>%
  dplyr::select(
    sample_id,
    Kids_First_Participant_ID,
    Kids_First_Biospecimen_ID
  )

#### Filter expression data ----------------------------------------------------
# Read in the stranded expression data file
expression_data <- readr::read_rds(file.path(
  root_dir,
  "data",
  "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
))

# Filter to Chordoma samples only -- we can use chordoma_df because it is subset to
# RNA-seq samples
expression_data <- expression_data %>%
  dplyr::select(intersect(
    chordoma_metadata$Kids_First_Biospecimen_ID,
    colnames(expression_data)
  )) %>%
  readr::write_rds(file.path(
    results_dir,
    "chordoma-only-gene-expression-rsem-fpkm-collapsed.stranded.rds"
  ))

#### Filter focal CN data ------------------------------------------------------

# Filter focal CN to Chordoma samples only
cn_metadata <- data.table::fread(file.path(
  root_dir,
  "data",
  "consensus_seg_annotated_cn_autosomes.tsv.gz"
)) %>%
  dplyr::left_join(chordoma_metadata,
    by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")
  ) %>%
  dplyr::select(
    gene_symbol,
    sample_id,
    Kids_First_Participant_ID,
    biospecimen_id,
    copy_number,
    ploidy,
    status
  ) %>%
  dplyr::filter(sample_id %in% chordoma_metadata$sample_id) %>%
  # Write to file
  readr::write_tsv(file.path(results_dir, "chordoma-only_cn_autosomes.tsv.gz"))
