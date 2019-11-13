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
  file.path(root_dir, "analyses", "focal-cn-file-preparation", "results")

plots_dir <-
  file.path(root_dir, "analyses", "focal-cn-file-preparation", "plots")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in data from tsv file (produced in `01-prepare-cn-file.R`)
cn_df <-
  readr::read_tsv(file.path(results_dir, "controlfreec_annotated_cn_autosomes.tsv.gz"))

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

#### Join data -----------------------------------------------------------------

# Add metadata to focal CN data.frame
cn_df_loss_metadata <- cn_df %>%
  dplyr::filter(status %in% c("loss")) %>%
  dplyr::inner_join(metadata,
                    by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(
    gene_symbol,
    biospecimen_id,
    status,
    copy_number,
    tumor_ploidy,
    Kids_First_Participant_ID,
    tumor_descriptor
  )

cn_df_non_loss_metadata <- cn_df %>%
  dplyr::filter(!(status %in% c("loss"))) %>%
  dplyr::inner_join(metadata,
                    by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(
    gene_symbol,
    biospecimen_id,
    status,
    copy_number,
    tumor_ploidy,
    Kids_First_Participant_ID,
    tumor_descriptor
  )

# Merge Focal CN data.frame with RNA expression data.frame using
# `merge_expression` custom function
rsem_combined_polyA_loss_df <-
  merge_expression(cn_df_loss_metadata,
                   rsem_expression_polyA,
                   metadata,
                   "annotated_cn_polyA_loss.tsv.gz")

rsem_combined_stranded_loss_df <-
  merge_expression(cn_df_loss_metadata,
                   rsem_expression_stranded,
                   metadata,
                   "annotated_cn_stranded_loss.tsv.gz")

rsem_combined_polyA_non_loss_df <-
  merge_expression(cn_df_non_loss_metadata,
                   rsem_expression_polyA,
                   metadata,
                   "annotated_cn_polyA_non_loss.tsv.gz")

rsem_combined_stranded_non_loss_df <-
  merge_expression(cn_df_non_loss_metadata,
                   rsem_expression_stranded,
                   metadata,
                   "annotated_cn_stranded_non_loss.tsv.gz")
#### Plot and Save -------------------------------------------------------------

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

# Plot and save using `plot_stacked_expression` custom function on polyA
# expression data
polyA_loss_plot_df <-
  plot_stacked_expression(
    rsem_combined_polyA_loss_df,
    "cn_loss_expression_polyA.png",
    "cn_loss_expression_per_gene_polyA.png"
  )

polyA_non_loss_plot_df <-
  plot_stacked_expression(
    rsem_combined_polyA_non_loss_df,
    "cn_non_loss_expression_polyA.png",
    "cn_non_loss_expression_per_gene_polyA.png"
  )

# Plot and save using `plot_stacked_expression` custom function on stranded
# expression data
stranded_loss_plot_df <-
  plot_stacked_expression(
    rsem_combined_stranded_loss_df,
    "cn_loss_expression_stranded.png",
    "cn_loss_expression_per_gene_stranded.png"
  )

stranded_non_loss_plot_df <-
  plot_stacked_expression(
    rsem_combined_stranded_non_loss_df,
    "cn_non_loss_expression_stranded.png",
    "cn_non_loss_expression_per_gene_stranded.png"
  )

# Plot and save scatterplot showing mean expression of deletions compared to 
# mean expression of non-deletions across genes 

mean_polyA_plot <-
  plot_mean_expression(
    polyA_loss_plot_df,
    polyA_non_loss_plot_df,
    "mean_scatterplot_polyA.png"
  )

mean_stranded_plot <-
  plot_mean_expression(
    stranded_loss_plot_df,
    stranded_non_loss_plot_df,
    "mean_scatterplot_stranded.png"
  )
