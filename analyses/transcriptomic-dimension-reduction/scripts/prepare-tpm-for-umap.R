###############################################################################################################
# SJ Spielman for CCDL, 2023
#
# This script collapses and prepares TPM data for UMAPs, where we run UMAP on two versions of the data:
# - One with all genes 
# - One without any of the mitochondrial genes. We do a quick re-scaling of TPM values after filtering out mito genes (`TPM = TPM/(upated TPM sum/1e6)`).
# To collapse genes, we read in the table exported by the `collapse-rnaseq` module (script linked below):
#    https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/059047a19025f4afd9143fceff174e6da3392747/analyses/collapse-rnaseq/01-summarize_matrices.R
# This table indicates which ensembl IDs were kept, so we'll filter to those. 
###############################################################################################################


# Directories and files -------

`%>%` <- dplyr::`%>%`
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

tpm_file <- file.path(root_dir, 
                      "data", 
                      "pbta-gene-expression-rsem-tpm.stranded.rds")
collapsed_file <- file.path(root_dir, 
                            "analyses", 
                            "collapse-rnaseq", 
                            "results", 
                            "pbta-gene-expression-rsem-fpkm-collapsed_table.stranded.rds")

output_dir <- file.path(root_dir,
                        "scratch", 
                        "transcriptomic-dimension-reduction")
if (!(dir.exists(output_dir))) {
  dir.create(output_dir, recursive = TRUE)
}
output_file_all <- file.path(output_dir,
                             "tpm-collapsed-all-stranded.rds")
output_file_nomito <- file.path(output_dir,
                                "tpm-collapsed-nomito-stranded.rds")

# Read input ------------

# Read TPM and separate ensembl and gene_symbol columns
tpm_df <- readr::read_rds(tpm_file) %>%
  tidyr::separate(gene_id, 
                  into = c("ensembl", "gene_symbol"),
                  sep = "_", 
                  extra = "merge") 

# Read in the collapsed df
collapsed_df <- readr::read_rds(collapsed_file)

# Filter to the collapsed gene set ------------
keep_ensembl <- collapsed_df$gene_id[collapsed_df$keep == "Yes"]
tpm_df_collapsed <- tpm_df %>%
  dplyr::filter(ensembl %in% keep_ensembl) %>%
  # no longer needed
  dplyr::select(-ensembl) 


# Create a version without mitochondrial genes -----------
# Filter out mitochondria genes and re-scale
tpm_df_collapsed_nomito <- tpm_df_collapsed %>%
  # make it long to facilitate filtering and renormalization
  tidyr::gather(
    key = "bs_ids",
    value = "tpm",
    tidyselect::starts_with("BS_")
  ) %>%
  # Remove mito
  dplyr::filter(!(stringr::str_starts(gene_symbol, "MT-"))) %>%
  # sample-level sum
  dplyr::group_by(bs_ids) %>%
  dplyr::mutate(total_tpm = (sum(tpm))/1e6) %>%
  # divde back out
  dplyr::mutate(tpm = tpm/(total_tpm)) %>%
  # remove total
  dplyr::select(-total_tpm) %>%
  # back to wide
  tidyr::spread(
    "bs_ids", "tpm"
  ) %>%
  # move gene_symbol into rownames for matrix conversion next
  tibble::column_to_rownames(var = "gene_symbol") 
  
# Reformat data frames back into matrices that are formatted as the transcriptomic-dimension-reduction module expects
# It expects gene symbols sownames, and BS ids solumn names, with TPM values.

tpm_collapsed_mat <- as.matrix(
  # This df still needs its rownames moved out
  tibble::column_to_rownames(tpm_df_collapsed, 
                             var = "gene_symbol") 
)
tpm_collapsed_nomito_mat <- as.matrix(tpm_df_collapsed_nomito)


# Export matrices to scratch
readr::write_rds(tpm_collapsed_mat, output_file_all)
readr::write_rds(tpm_collapsed_nomito_mat, output_file_nomito)




