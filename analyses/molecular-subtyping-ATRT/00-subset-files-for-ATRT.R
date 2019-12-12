# This script subsets the focal copy number, RNA expression, tumor mutation 
# burden and histologies` files to include only ATRT samples.

# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/molecular-subtyping-ATRT/00-subset-files-for-ATRT.R'


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
  file.path(root_dir, "analyses", "molecular-subtyping-ATRT", "atrt-subset")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv"))

# Select wanted columns in metadata for merging and assign to a new object
select_metadata <- metadata %>%
  dplyr::select(sample_id,
                Kids_First_Participant_ID,
                Kids_First_Biospecimen_ID)

# Read in ssGSEA pathway information
ssGSEA <-
  as.data.frame(readr::read_rds(
    file.path(
      root_dir,
      "analyses",
      "ssgsea-hallmark",
      "results",
      "GeneSetExpressionMatrix.RDS"
    )
  ))

# Read in RNA expression data
stranded_expression <-
  readr::read_rds(
    file.path(
      root_dir,
      "data",
      "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
    )
  )

# Read in focal CN data
## TODO: This section will be updated to read in focal CN data derived from
##       copy number consensus calls.
cn_df <- readr::read_tsv(
  file.path(
    root_dir,
    "analyses",
    "focal-cn-file-preparation",
    "results",
    "controlfreec_annotated_cn_autosomes.tsv.gz"
  )
)

# Read in consensus mutation data
tmb_df <-
  data.table::fread(file.path(root_dir,
                              "data",
                              "pbta-snv-consensus-mutation-tmb.tsv"))

## TODO: Read in the SV data/GISTIC output to evaluate the chr22q loss 
#        associated with SMARB1 deletions

#### Filter metadata -----------------------------------------------------------

# Filter metadata for `ATRT` and define `location_summary` based on values in
# `primary_site`
atrt_df <- metadata %>%
  dplyr::filter(short_histology == "ATRT",
                experimental_strategy == "RNA-Seq")

# Write to file
readr::write_tsv(atrt_df, file.path(results_dir, "atrt_histologies.tsv"))

#### Filter expression data ----------------------------------------------------

# Filter to ATRT samples only -- we can use atrt_df because it is subset to
# RNA-seq samples
stranded_expression <- stranded_expression %>%
  dplyr::select(atrt_df$Kids_First_Biospecimen_ID)

# Log2 transformation
norm_expression <- log2(stranded_expression + 1)

# normData mean and sd
norm_expression_means <- rowMeans(norm_expression, na.rm = TRUE)
norm_expression_sd <- apply(norm_expression, 1, sd, na.rm = TRUE)

# Subtract mean
expression_zscored <-
  sweep(norm_expression, 1, norm_expression_means, FUN = "-")

# Divide by SD remove NAs and Inf values from zscore for genes with 0 in normData
expression_zscored <-
  sweep(expression_zscored, 1, norm_expression_sd, FUN = "/") %>%
  dplyr::na_if(Inf) %>%
  na.omit()

# Save matrix with all genes to file for downstream plotting
readr::write_rds(expression_zscored, 
                 file.path(results_dir, "atrt_zscored_expression.RDS"))

#### Filter focal CN data ------------------------------------------------------

# Filter focal CN to ATRT samples only 
cn_metadata <- cn_df %>%
  dplyr::left_join(select_metadata,
                   by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(gene_symbol,
                sample_id,
                Kids_First_Participant_ID,
                biospecimen_id,
                status) %>%
  dplyr::filter(sample_id %in% atrt_df$sample_id) %>%
  dplyr::distinct()

# Write to file
readr::write_tsv(cn_metadata, file.path(results_dir, "atrt_focal_cn.tsv.gz"))

#### Filter ssGSEA data --------------------------------------------------------

# Transpose
transposed_ssGSEA <- t(ssGSEA) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID")

# Filter for `sample_id` values found in the metadata filtered for ATRT samples
transposed_ssGSEA <- transposed_ssGSEA %>%
  dplyr::left_join(select_metadata, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::filter(sample_id %in% atrt_df$sample_id) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID") %>%
  dplyr::select(-c("sample_id", "Kids_First_Participant_ID"))

# Write to file
readr::write_rds(transposed_ssGSEA, file.path(results_dir, "atrt_ssgsea.RDS"))

#### Filter tumor mutation burden data -----------------------------------------

tmb_df <- tmb_df %>%
  dplyr::select(Tumor_Sample_Barcode, tmb) %>%
  dplyr::left_join(select_metadata,
                    by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::filter(sample_id %in% atrt_df$sample_id)

# Write to file 
readr::write_tsv(tmb_df, file.path(results_dir, "atrt_tmb.tsv"))
