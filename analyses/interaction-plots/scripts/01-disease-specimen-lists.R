# Create sample lists by disease type
#
# JA Shapiro for ALSF - CCDL
#
# 2019
#
# Generates a table of sample lists for disease type analysis

# Option descriptions
#
# --metadata : Relative file path to metadata with sample information.
#     File path is given from top directory of 'OpenPBTA-analysis'.
# --specimen_list: A file of specimens to include. Ideally, this list should consist
#     of independent samples (at most one from each individual). File path is given
#     from the top directory of 'OpenPBTA-analysis'.
# --disease: The name of the disease to be analyzed (in quotes).
#     "All" for any disease type.
# --outfile: The file path for output (relative to top of OpenPBTA-analysis)
#
# Command line example:
#
# Rscript analyses/interaction-plots/01-process-mutations.R \
#   --metadata data/pbta-histologies.tsv
#   --specimen_list data/independent-specimens.wgs.primary.tsv
#   --disease "Medulloblastoma"
#   --outfile scratch/medulloblastoma_specimens.tsv

#### Initial Set Up
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Load libraries:
library(optparse)


#### Set up options
option_list <- list(
  make_option(
    opt_str = "--metadata",
    type = "character",
    default = file.path("data", "pbta-histologies.tsv"),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--specimen_list",
    type = "character",
    default = file.path(
      "data",
      "independent-specimens.wgs.primary.tsv"
    ),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            to MAF file to be analyzed. Can be .gz compressed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--outfile",
    type = "character",
    default = file.path("analyses", "interaction-plots", "results", "cooccurence.tsv"),
    help = "Relative file path (from top directory of 'OpenPBTA-analysis')
            where output table will be placed.",
    metavar = "character"
  ),
  make_option(
    opt_str = "--disease",
    type = "character",
    default = "All",
    help = "Disease type, as found in the `disease_type_new` column of the histology file.",
    metavar = "character"
  )
)

# Parse options
opts <- parse_args(OptionParser(option_list = option_list))

# file locations
meta_file <- file.path(root_dir, opts$metadata)
specimen_file <- file.path(root_dir, opts$specimen_list)
out_file <- file.path(root_dir, opts$outfile)


#### Read files

meta_df <- readr::read_tsv(meta_file, col_types = readr::cols(), guess_max = 10000)
specimen_df <- readr::read_tsv(specimen_file, col_types = readr::cols())


### reduce metadata to specimen list and disease type

disease_df <- meta_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% specimen_df$Kids_First_Biospecimen_ID) %>%
  dplyr::filter(tolower(opts$disease) == "all" |
    tolower(disease_type_new) == tolower(opts$disease)) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)

readr::write_tsv(disease_df, out_file)
