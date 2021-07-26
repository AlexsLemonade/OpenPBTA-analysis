# This script examines RNA-seq expression levels of genes that are called as
# deletions.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript --vanilla rna-expression-validation.R \
#   --annotated_cnv_file analyses/focal-cn-file-preparation/results/cnvkit_annotated_cn_autosomes.tsv.gz \
#   --expression_file data/gene-expression-rsem-tpm-collapsed.rds \
#   --independent_specimens_file data/independent-specimens.wgswxs.primary.tsv \
#   --metadata  data/histologies.tsv \
#   --goi_list analyses/oncoprint-landscape/driver-lists/brain-goi-list-long.txt \
#   --filename_lead "cnvkit_annotated_cn_autosomes"

#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Command line options ------------------------------------------------------

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("--annotated_cnv_file"),
    type = "character",
    default = NULL,
    help = "file path to file that contains annoatated CNV information"
  ),
  optparse::make_option(
    c("--expression_file"),
    type = "character",
    default = NULL,
    help = "file path to RDS file that contains gene expression information"
  ),
  optparse::make_option(
    c("--independent_specimens_file"),
    type = "character",
    default = NULL,
    help = "file path to tsv file that contains list of independent specimen ids"
  ),
  optparse::make_option(
    c("--metadata"),
    type = "character",
    default = NULL,
    help = "file path to histologies.tsv"
  ),
  optparse::make_option(
    c("--goi_list"),
    type = "character",
    default = NULL,
    help = "file path to tsv file that contains list of driver genes"
  ),
  optparse::make_option(
    c("--filename_lead"),
    type = "character",
    default = "annotated_expression",
    help = "used in file names"
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to results and plots directory
results_dir <- file.path(root_dir, "scratch")

plots_dir <-
  file.path(root_dir, "analyses", "focal-cn-file-preparation", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in data from tsv file (produced in `03-prepare-cn-file.R`)
cn_df <- readr::read_tsv(opt$annotated_cnv_file)

# Read in RNA-seq expression data
expression_matrix <- readr::read_rds(opt$expression_file)

# Read in metadata
metadata <- readr::read_tsv(opt$metadata,
                            col_types = readr::cols(molecular_subtype = readr::col_character()))

#### Custom Functions ----------------------------------------------------------

# Source script with custom functions
source(
  file.path(
    root_dir,
    "analyses",
    "focal-cn-file-preparation",
    "util",
    "rna-expression-functions.R"
  )
)

#### Get rid of ambiguous and non-tumor samples --------------------------------

# Below code is adapted from: analyses/oncoprint-landscape/00-map-to-sample_id.R
# An ambiguous sample_id will have more than 2 rows associated with it in the
# histologies file when looking at tumor samples -- that means we won't be able
# to determine when an WGS/WXS assay maps to an RNA-seq assay for the purpose of
# plotting
ambiguous_sample_ids <- metadata %>%
  dplyr::filter(
    sample_type == "Tumor",
    composition == "Solid Tissue" | composition == "Bone Marrow"
  ) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::tally() %>%
  dplyr::filter(n > 2) %>%
  dplyr::pull(sample_id)

ambiguous_biospecimens <- metadata %>%
  dplyr::filter(sample_id %in% ambiguous_sample_ids) %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

# We're only going to look at tumor samples
not_tumor_biospecimens <- metadata %>%
  dplyr::filter(
    sample_type != "Tumor",
    composition != "Solid Tissue" & composition != "Bone Marrow"
  ) %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

biospecimens_to_remove <- unique(c(
  ambiguous_biospecimens,
  not_tumor_biospecimens
))

# Filter the CN data
cn_df <- cn_df %>%
  dplyr::filter(!(biospecimen_id %in% biospecimens_to_remove)) %>%
  dplyr::select(biospecimen_id, status, copy_number, gene_symbol)

#### Determine independent specimens -------------------------------------------

# Below code is adapted from: analyses/oncoprint-landscape/00-map-to-sample_id.R
# Read in the primary indepedent specimens file
ind_biospecimen <-
  readr::read_tsv(opt$independent_specimens_file) %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

# Filter the CN data, to only include biospecimen identifiers in the
# independent file
cn_df <- cn_df %>%
  dplyr::filter(biospecimen_id %in% ind_biospecimen) %>%
  dplyr::distinct() %>%
  dplyr::group_by(biospecimen_id, gene_symbol) %>%
  # Collapse status and copy number -- anything with conflicting (e.g., more 
  # than one) values will contain a comma in status and/or copy number
  dplyr::summarize(status = paste(sort(unique(status)),
                                  collapse = ", "),
                   copy_number = paste(sort(unique(copy_number)),
                                       collapse = ", "))

# TODO: write cases where there is conflicting evidence re: status and copy
# number to a separate file -- this would be the same regardless of the 
# expression data being used 
# ambiguous_cn_status_df <- cn_df %>%
#   dplyr::filter(stringr::str_detect(status, ",") |
#                   stringr::str_detect(copy_number, ","))

# This is removing any instances where the status or copy number do not agree
cn_df <-cn_df %>%
  dplyr::filter(stringr::str_detect(status, ",", negate = TRUE) |
                  stringr::str_detect(copy_number, ",", negate = TRUE))

# For the RNA-seq samples, we need to map from the sample identifier
# associated with the independent specimen and back to a biospecimen ID
ind_sample_id <- metadata %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% ind_biospecimen) %>%
  dplyr::pull(sample_id)

# Get the corresponding biospecimen ID to be used
rnaseq_ind <- metadata %>%
  dplyr::filter(
    sample_id %in% ind_sample_id,
    experimental_strategy == "RNA-Seq"
  ) %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

#### Filter and Join data ------------------------------------------------------

# Filter expression data to include independent sample and calculate z-scores by
# gene using `calculate_z_score` custom function
expression_zscore_df <- calculate_z_score(expression_matrix, rnaseq_ind)

# Merge Focal CN data.frame with RNA expression data.frame using
# `merge_expression` custom function
expression_cn_combined_df <-
  merge_expression(
    cn_df,
    expression_zscore_df,
    metadata,
    opt$filename_lead
  )

# Read in gene list
goi_list <-
  read.delim(opt$goi_list,
    sep = "\t",
    header = FALSE,
    as.is = TRUE
  )

# Add `is_driver_gene` column denoting whether or not the gene is found in the
# `goi_list` of driver genes
filtered_expression_cn_combined_df <- expression_cn_combined_df %>%
  dplyr::mutate(is_driver_gene = ifelse(gene_id %in% goi_list$V1, "Yes", "No"))

#### Plot and Save -------------------------------------------------------------

# Filter for loss CN calls
cn_loss_df <- filtered_expression_cn_combined_df %>%
  dplyr::filter(status == "loss")

# Filter for neutral CN calls
cn_neutral_df <- filtered_expression_cn_combined_df %>%
  dplyr::filter(is.na(status))

# Filter for copy number equal to 0
cn_zero_df <- filtered_expression_cn_combined_df %>%
  dplyr::filter(copy_number == 0)

# Plot and save using `plot_stacked_expression` custom function on expression
# data
cn_loss_plot <-
  plot_stacked_expression(
    cn_loss_df,
    cn_neutral_df,
    cn_zero_df,
    opt$filename_lead
  )

# Plot and save scatterplot showing mean expression of deletions compared to
# mean expression of non-deletions across genes
cn_mean_plot <-
  plot_mean_expression(
    cn_loss_df,
    cn_neutral_df,
    cn_zero_df,
    opt$filename_lead
  )

