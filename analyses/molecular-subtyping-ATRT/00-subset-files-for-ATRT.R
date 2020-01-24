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

# Read in GSVA pathway scores
gsva_scores <-
  as.data.frame(readr::read_tsv(
    file.path(
      root_dir,
      "analyses",
      "gene-set-enrichment-analysis",
      "results",
      "gsva_scores_stranded.tsv"
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
    "cnvkit_annotated_cn_autosomes.tsv.bz2"
  )
)

# Read in consensus mutation data
tmb_df <-
  data.table::fread(file.path(root_dir,
                              "data",
                              "pbta-snv-consensus-mutation-tmb-all.tsv"))

# Read in GISTIC `broad_values_by_arm.txt` file
gistic_df <-
  data.table::fread(unzip(
    file.path(root_dir, "data", "pbta-cnv-cnvkit-gistic.zip"),
    files = file.path(
      "2019-12-10-gistic-results-cnvkit",
      "broad_values_by_arm.txt"
    ),
    exdir = file.path(root_dir, "scratch")
  ), data.table = FALSE)

#### Filter metadata -----------------------------------------------------------

atrt_df <- metadata %>%
  dplyr::filter(short_histology == "ATRT",
                sample_type == "Tumor",
                composition == "Solid Tissue")

#### Filter expression data ----------------------------------------------------

# Filter to ATRT samples only -- we can use atrt_df because it is subset to
# RNA-seq samples
stranded_expression <- stranded_expression %>%
  dplyr::select(intersect(atrt_df$Kids_First_Biospecimen_ID,
                          colnames(stranded_expression)))

# Log2 transformation
norm_expression <- log2(stranded_expression + 1)

# Save matrix with all genes to file for downstream plotting
readr::write_rds(norm_expression,
                 file.path(results_dir, "atrt_log_expression.RDS"))

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

#### Filter GSVA data ----------------------------------------------------------

# Filter for `sample_id` values found in the metadata filtered for ATRT samples
filtered_gsva_scores <- gsva_scores %>%
  dplyr::left_join(select_metadata, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::filter(sample_id %in% atrt_df$sample_id) %>%
  dplyr::select(-Kids_First_Participant_ID)

# Write to file
readr::write_tsv(filtered_gsva_scores, file.path(results_dir, "atrt_gsva.tsv"))

#### Filter tumor mutation burden data -----------------------------------------

tmb_df <- tmb_df %>%
  dplyr::select(Tumor_Sample_Barcode, tmb) %>%
  dplyr::left_join(select_metadata,
                   by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::filter(sample_id %in% atrt_df$sample_id)

# Write to file
readr::write_tsv(tmb_df, file.path(results_dir, "atrt_tmb.tsv"))

#### Filter GISTIC data --------------------------------------------------------

gistic_df <- gistic_df %>%
  tibble::column_to_rownames("Chromosome Arm") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% atrt_df$Kids_First_Biospecimen_ID) %>%
  dplyr::left_join(select_metadata, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(sample_id,
                Kids_First_Biospecimen_ID,
                `22q`) # Select only the chromosome arm we are interested in

# Write to file
readr::write_tsv(gistic_df,
                 file.path(results_dir, "atrt_gistic_broad_values.tsv"))
