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
  function (copy_number_df, expression_df, metadata, filename){
    # Given the focal copy number data.frame already annotated with the
    # metadata, the RNA-seq expression data.frame, and the metadata, combine
    # the data into one data.frame and save data.frame as tsv file. 
    #
    # Args:
    #   copy_number_df: focal copy number data.frame annotated with the metadata
    #   expression_df: RNA-seq expression data.frame
    #   metadata: the relevant metadata data.frame
    #   filename: filename of the output tsv file
    #
    # Returns:
    #   combined_subset_df: data.frame with information from the focal CN, the
    #                       RNA-seq expression, and metadata data.frames, which
    #                       has been subsetted for plotting

    # Change expression data.frame to long tidy format
    long_expression <- expression_df %>%
      tidyr::gather(biospecimen_id, expression_value, -gene_id) %>%
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
      ) %>%
      dplyr::rename(biospecimen_id_cn = biospecimen_id.x, 
                    biospecimen_id_expression = biospecimen_id.y) %>%
      dplyr::mutate(log_expression_value = (log2(expression_value) + 1)) %>%
      dplyr::mutate(expression_class = dplyr::case_when(
        expression_value == 0 ~ "0",
        expression_value > 0 & expression_value <= 10 ~ "0-10",
        expression_value > 10 & expression_value <= 100 ~ "10-100",
        expression_value > 100 ~ ">100"))

    # Save results
    readr::write_tsv(combined_df, file.path(results_dir, filename))

    # Subset data.frame for plotting
    combined_subset_df <- combined_df %>%
      dplyr::filter(!(expression_class == ">100"))
  }

#### Join data -----------------------------------------------------------------

# Add metadata to focal CN data.frame
cn_df_metadata <- cn_df %>%
  dplyr::filter(label %in% c("Hom_Deletion", "Hem_Deletion")) %>%
  dplyr::inner_join(metadata,
                    by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(
    gene_symbol,
    biospecimen_id,
    label,
    copy_number,
    tumor_ploidy,
    Kids_First_Participant_ID,
    tumor_descriptor
  )

# Merge Focal CN data.frame with RNA expression data.frame using
# `merge_expression` custom function
rsem_combined_polyA_df <-
  merge_expression(cn_df_metadata,
                   rsem_expression_polyA,
                   metadata,
                   "annotated_cn_polyA.tsv.gz")

rsem_combined_stranded_df <-
  merge_expression(cn_df_metadata,
                   rsem_expression_stranded,
                   metadata,
                   "annotated_cn_stranded.tsv.gz")

#### Plot and Save -------------------------------------------------------------

png(
  file.path(plots_dir, "cn_loss_expression_polyA.png"),
  width = 1500,
  height = 900,
  res = 120
)
ggplot2::ggplot(rsem_combined_polyA_df,
                ggplot2::aes(x = label, fill = expression_class)) +
  ggplot2::geom_bar()
dev.off()

png(
  file.path(plots_dir, "cn_loss_expression_stranded.png"),
  width = 1500,
  height = 900,
  res = 120
)
ggplot2::ggplot(rsem_combined_stranded_df,
                ggplot2::aes(x = label, fill = expression_class)) +
  ggplot2::geom_bar()
dev.off()
