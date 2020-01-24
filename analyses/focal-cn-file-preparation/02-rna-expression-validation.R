# This script examines RNA-seq expression levels of genes that are called as
# deletions.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/focal-cn-file-preparation/rna-expression-validation.R'

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
  file.path(root_dir, "scratch")

plots_dir <-
  file.path(root_dir, "analyses", "focal-cn-file-preparation", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in data from tsv file (produced in `01-prepare-cn-file.R`)
cn_df <-
  readr::read_tsv(file.path(root_dir,
                            "analyses",
                            "focal-cn-file-preparation",
                            "results",
                            "controlfreec_annotated_cn_autosomes.tsv.bz2"))

# Read in RNA-seq expression data
rsem_expression_polyA <-
  readr::read_rds(file.path(root_dir,
                            "data",
                            "pbta-gene-expression-rsem-fpkm.polya.rds"))

rsem_expression_stranded <-
  readr::read_rds(file.path(
    root_dir,
    "data",
    "pbta-gene-expression-rsem-fpkm.stranded.rds"
  ))

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv"))

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
  dplyr::filter(sample_type == "Tumor",
         composition == "Solid Tissue") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::tally() %>%
  dplyr::filter(n > 2) %>%
  dplyr::pull(sample_id)

ambiguous_biospecimens <- metadata %>%
  dplyr::filter(sample_id %in% ambiguous_sample_ids) %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

# We're only going to look at tumor samples
not_tumor_biospecimens <- metadata %>%
  dplyr::filter(sample_type != "Tumor",
         composition != "Solid Tissue") %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

biospecimens_to_remove <- unique(c(ambiguous_biospecimens,
                                   not_tumor_biospecimens))

# Filter the CN data 
cn_df <- cn_df %>%
  dplyr::filter(!(biospecimen_id %in% biospecimens_to_remove))


#### Determine independent specimens -------------------------------------------

# Below code is adapted from: analyses/oncoprint-landscape/00-map-to-sample_id.R
# Read in the primary indepedent specimens file
ind_biospecimen <-
  readr::read_tsv(file.path(root_dir,
                            "data",
                            "independent-specimens.wgswxs.primary.tsv")) %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

# Filter the CN data, to only include biospecimen identifiers in the 
# independent file
cn_df <- cn_df %>%
  dplyr::filter(biospecimen_id %in% ind_biospecimen)
  
# For the RNA-seq samples, we need to map from the sample identifier
# associated with the independent specimen and back to a biospecimen ID
ind_sample_id <- metadata %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% ind_biospecimen) %>%
  dplyr::pull(sample_id)
  
# Get the corresponding biospecimen ID to be used
rnaseq_ind <- metadata %>%
  dplyr::filter(sample_id %in% ind_sample_id,
                experimental_strategy == "RNA-Seq") %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

#### Filter and Join data ------------------------------------------------------

# Filter expression data to include independent sample and calculate z-scores by 
# gene using `calculate_z_score` custom function
expression_df_polyA <- calculate_z_score(rsem_expression_polyA, rnaseq_ind)
expression_df_stranded <- calculate_z_score(rsem_expression_stranded,
                                            rnaseq_ind)

# Merge Focal CN data.frame with RNA expression data.frame using
# `merge_expression` custom function
rsem_combined_polyA_df <-
  merge_expression(cn_df,
                   expression_df_polyA,
                   metadata,
                   "annotated_focal_cn_expression_polyA.tsv.gz")

rsem_combined_stranded_df <-
  merge_expression(cn_df,
                   expression_df_stranded,
                   metadata,
                   "annotated_focal_cn_expression_stranded.tsv.gz")

# Read in gene list
goi_list <-
  read.delim(
    file.path(
      root_dir,
      "analyses",
      "oncoprint-landscape",
      "driver-lists",
      "brain-goi-list-long.txt"
    ),
    sep = "\t",
    header = FALSE,
    as.is = TRUE
  )

# Add `is_driver_gene` column denoting whether or not the gene is found in the
# `goi_list` of driver genes
rsem_combined_polyA_df <- rsem_combined_polyA_df %>%
  dplyr::mutate(is_driver_gene = ifelse(gene_id %in% goi_list$V1, "Yes", "No"))

rsem_combined_stranded_df <- rsem_combined_stranded_df %>%
  dplyr::mutate(is_driver_gene = ifelse(gene_id %in% goi_list$V1, "Yes", "No"))

#### Plot and Save -------------------------------------------------------------

# Filter for loss CN calls
cn_loss_polyA_df <- rsem_combined_polyA_df %>%
  dplyr::filter(status == "loss")

cn_loss_stranded_df <- rsem_combined_stranded_df %>%
  dplyr::filter(status == "loss")

# Filter for neutral CN calls
cn_neutral_polyA_df <- rsem_combined_polyA_df %>%
  dplyr::filter(is.na(status))

cn_neutral_stranded_df <- rsem_combined_stranded_df %>%
  dplyr::filter(is.na(status))

# Filter for copy number equal to 0
cn_zero_polyA_df <- rsem_combined_polyA_df %>%
  dplyr::filter(copy_number == 0)

cn_zero_stranded_df <- rsem_combined_stranded_df %>%
  dplyr::filter(copy_number == 0)

# Plot and save using `plot_stacked_expression` custom function on expression
# data
polyA_loss_plot_df <-
  plot_stacked_expression(
    cn_loss_polyA_df,
    cn_neutral_polyA_df,
    cn_zero_polyA_df,
    "focal_cn_stacked_expression_polyA.png"
  )

stranded_loss_plot_df <-
  plot_stacked_expression(
    cn_loss_stranded_df,
    cn_neutral_stranded_df,
    cn_zero_stranded_df,
    "focal_cn_stacked_expression_stranded.png"
  )

# Plot and save scatterplot showing mean expression of deletions compared to 
# mean expression of non-deletions across genes 

mean_polyA_plot <-
  plot_mean_expression(
    cn_loss_polyA_df,
    cn_neutral_polyA_df,
    cn_zero_polyA_df,
    "loss_correlation_scatterplot_polyA.png",
    "zero_correlation_scatterplot_polyA.png"
  )

mean_stranded_plot <-
  plot_mean_expression(
    cn_loss_stranded_df,
    cn_neutral_stranded_df,
    cn_zero_stranded_df,
    "loss_correlation_scatterplot_stranded.png",
    "zero_correlation_scatterplot_stranded.png"
  )
