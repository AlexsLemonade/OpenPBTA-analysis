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

output_dir <- file.path(root_dir,
                        "scratch", 
                        "transcriptomic-dimension-reduction")
if (!(dir.exists(output_dir))) {
  dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir,
                         "tpm-nomito-stranded.rds")

# Remove mitogenes from TPM -------


# Read in and filter out mito genes
tpm_df <- readr::read_rds(tpm_file) %>%
  # separate ensembl and gene_symbol columns so we can filter out mito genes
  tidyr::separate(gene_id, 
                  into = c("ensembl", "gene_symbol"),
                  sep = "_", 
                  extra = "merge", 
                  # we'll want to restore `gene_id`
                  remove = FALSE) %>%
  # remove mito genes
  dplyr::filter(!(stringr::str_starts(gene_symbol, "MT-"))) %>%
  # remove extra separated columns now that we have filtered
  dplyr::select(-ensembl, -gene_symbol) %>%
  # make the df long to facilitate renormalization
  tidyr::gather(
      key = "bs_ids",
      value = "tpm",
      tidyselect::starts_with("BS_")
      ) 

# Re-normalize now that mito genes are removed
tpm_df <- tpm_df %>%
  # calculate sample-level sum
  dplyr::group_by(bs_ids) %>%
  dplyr::mutate(total_tpm = (sum(tpm))/1e6) %>%
  # divde back out for new tpm
  dplyr::mutate(tpm = tpm/(total_tpm)) %>%
  # remove total
  dplyr::select(-total_tpm) %>%
  # back to wide, as expected by the `transcriptomic-dimension-reduction` module
  tidyr::spread(
    "bs_ids", "tpm"
  ) 
  
# Export matrix to scratch
readr::write_rds(tpm_df, output_file)




