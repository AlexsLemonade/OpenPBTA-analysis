# Bethell and Taroni for CCDL 2019
# Subset the gene expression files (kallisto and RSEM) by the RNA library 
# strategy that is contained in the histologies TSV
#
# Usage:
# Rscript --vanilla 01-subset-expression-by-strategy.R

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# we consider these 
output_dir <- file.path(root_dir, "scratch")

#### Identify samples for each strategy ----------------------------------------

metadata_file <- file.path(root_dir, "data", "pbta-histologies.tsv")
metadata_df <- readr::read_tsv(metadata_file)

stranded_samples <- metadata_df %>%
  dplyr::filter(RNA_library == "stranded") %>%
  dplyr::pull(Kids_First_Biospecimen_ID)

polyA_samples <- metadata_df %>%
  dplyr::filter(RNA_library == "poly-A") %>%
  dplyr::pull(Kids_First_Biospecimen_ID)


#### RSEM ----------------------------------------------------------------------

rsem_all_file <- file.path(root_dir, "data", 
                           "pbta-gene-expression-rsem.fpkm.rds")
rsem <- readr::read_rds(rsem_all_file) 

# this handles errors related to 'missing' samples
rsem_stranded_samples <- intersect(colnames(rsem), stranded_samples)
rsem_polyA_samples <- intersect(colnames(rsem), polyA_samples)

rsem_stranded <- rsem %>%
  dplyr::select(gene_id, rsem_stranded_samples) %>%
  readr::write_rds(file.path(output_dir, 
                             "pbta-gene-expression-rsem_stranded.fpkm.rds"))

rsem_polyA <- rsem %>%
  dplyr::select(gene_id, rsem_polyA_samples) %>%
  readr::write_rds(file.path(output_dir, 
                             "pbta-gene-expression-rsem_polyA.fpkm.rds"))

#### kallisto ------------------------------------------------------------------

kallisto_all_file <- file.path(root_dir, "data", 
                               "pbta-gene-expression-kallisto.rds")
kallisto <- readr::read_rds(kallisto_all_file) 

# this handles errors related to 'missing' samples
kallisto_stranded_samples <- intersect(colnames(kallisto), stranded_samples)
kallisto_polyA_samples <- intersect(colnames(kallisto), polyA_samples)

kallisto_stranded <- kallisto %>%
  dplyr::select(transcript_id, gene_id, kallisto_stranded_samples) %>%
  readr::write_rds(file.path(output_dir, 
                             "pbta-gene-expression-kallisto_stranded.rds"))

kallisto_polyA <- kallisto %>%
  dplyr::select(transcript_id, gene_id, kallisto_polyA_samples) %>%
  readr::write_rds(file.path(output_dir, 
                             "pbta-gene-expression-kallisto_polyA.rds"))
