# This script examines RNA-seq expression levels of genes that are called as
# deletions.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/oncoprint-landscape/rna-expression-validation.R'

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
cn_df <- readr::read_tsv(file.path(results_dir, "annotated_cn.tsv"))

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

#### Custom Function -----------------------------------------------------------

merge_expression <-
  function(copy_number_df, expression_df, metadata) {
    # Given the focal copy number data.frame already annotated with the metadata,
    # the RNA-seq expression data.frame, and the metadata, combine the data into one
    # data.frame and plot expression values.
    #
    # Args:
    #   copy_number_df: focal copy number data.frame annotated with the metadata
    #   expression_df: RNA-seq expression data.frame
    #   metadata: the relevant metadata data.frame
    #
    # Returns:
    #   combined_df: data.frame with information from the focal CN, the
    #                RNA-seq expression, and metadata data.frames
    
    # Change expression data.frame to long tidy format
    long_expression <- expression_df %>%
      tidyr::gather(biospecimen_id, expression_value,-gene_id) %>%
      dplyr::distinct()
    long_expression$gene_id <-
      gsub(".*\\_", "", long_expression$gene_id)
    
    # Annotate expression data with metadata
    expression_metadata <- long_expression %>%
      dplyr::inner_join(metadata,
                        by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
      dplyr::select(
        biospecimen_id,
        expression_value,
        gene_id,
        Kids_First_Participant_ID,
        tumor_descriptor
      ) %>%
      dplyr::distinct()
    
    # Merge Focal CN data.frame with RNA expression data.frame
    combined_df <- copy_number_df %>%
      dplyr::distinct() %>%
      dplyr::inner_join(
        expression_metadata,
        by = c(
          "Kids_First_Participant_ID",
          "gene_symbol" = "gene_id",
          "tumor_descriptor"
        )
      )
  }

#### Join data -----------------------------------------------------------------

# Add metadata to focal CN data.frame
cn_df_metadata <- cn_df %>%
  dplyr::inner_join(metadata,
                    by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(
    gene_symbol,
    biospecimen_id,
    label,
    copy_number,
    Kids_First_Participant_ID,
    tumor_descriptor,
    tumor_ploidy
  ) %>%
  dplyr::filter(label %in% c("Hom_Deletion", "Hem_Deletion"))

# Merge Focal CN data.frame with RNA expression data.frame using
# `merge_expression` custom function
rsem_combined_polyA_df <-
  merge_expression(cn_df_metadata, rsem_expression_polyA, metadata)

rsem_combined_stranded_df <-
  merge_expression(cn_df_metadata, rsem_expression_stranded, metadata)

#### Plot and Save -------------------------------------------------------------

png(
  file.path(plots_dir, "cn_loss_expression_polyA.png"),
  width = 1200,
  height = 600
)
ggplot2::ggplot(rsem_combined_polyA_df,
                ggplot2::aes(x = gene_symbol, y = expression_value)) +
  ggplot2::geom_bar(stat = "identity", ggplot2::aes(col = label)) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      size = 0.2
    ),
    title = ggplot2::element_text("Focal CN polyA expression")
  )
dev.off()

png(
  file.path(plots_dir, "cn_loss_expression_stranded.png"),
  width = 1200,
  height = 600
)
ggplot2::ggplot(rsem_combined_stranded_df,
                ggplot2::aes(x = gene_symbol, y = expression_value)) +
  ggplot2::geom_point(ggplot2::aes(col = label), size = 0.6, alpha = 0.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      size = 0.5
    ),
    title = ggplot2::element_text("Focal CN stranded expression")
  )
dev.off()
