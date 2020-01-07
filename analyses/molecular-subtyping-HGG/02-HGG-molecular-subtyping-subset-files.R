# This script subsets the focal copy number, RNA expression, fusion and 
# histologies` files to include only High-grade glioma samples.

# Chante Bethell for CCDL 2020
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/molecular-subtyping-ATRT/02-HGG-molecular-subtyping-subset-files.R'


#### Set Up --------------------------------------------------------------------

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to subset directory
subset_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-HGG", "hgg-subset")

if (!dir.exists(subset_dir)) {
  dir.create(subset_dir)
}

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv"))

# Select wanted columns in metadata for merging and assign to a new object
select_metadata <- metadata %>%
  dplyr::select(sample_id,
                Kids_First_Participant_ID,
                Kids_First_Biospecimen_ID,
                glioma_brain_region,
                age_at_diagnosis_days)

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

# Read in fusion data
fusion_df <- readr::read_tsv(
  file.path(root_dir, "data", "pbta-fusion-putative-oncogenic.tsv"))

# Read in output file from `01-HGG-molecular-subtyping-defining-lesions.Rmd`
hgg_lesions_df <- readr::read_tsv(
  file.path(
    root_dir,
    "analyses",
    "molecular-subtyping-HGG",
    "results",
    "HGG_defining_lesions.tsv"
  )
)

#### Filter HGG defining lesions data.frame ------------------------------------

# Filter the output file from `01-HGG-molecular-subtyping-defining-lesions.Rmd`
# for samples classified and reclassified as HGG
hgg_lesions_df <- hgg_lesions_df %>%
  dplyr::filter(
    disease_type_new == "High-grade glioma" |
      grepl("High-grade glioma", disease_type_reclassified)
  )

#### Filter metadata -----------------------------------------------------------

# Filter metadata for `High-grade glioma` and samples that should be classified
# as High-grade glioma based on defining lesions and are stranded RNA-seq
# samples.
hgg_metadata_df <- metadata %>%
  dplyr::filter(
    disease_type_new == "High-grade glioma" |
      sample_id %in% hgg_lesions_df$sample_id,
    experimental_strategy == "RNA-Seq",
    RNA_library == "stranded"
  )

# Write to file
readr::write_tsv(hgg_metadata_df, file.path(subset_dir, "hgg_histologies.tsv"))

#### Filter expression data ----------------------------------------------------

# Filter to HGG samples only -- we can use hgg_metadata_df because it is subset
# to RNA-seq samples
stranded_expression <- stranded_expression %>%
  dplyr::select(hgg_metadata_df$Kids_First_Biospecimen_ID)

# Log2 transformation
norm_expression <- log2(stranded_expression + 1)

# Scale does column centering, so we transpose first
long_zscored_expression <- scale(t(norm_expression), 
                                  center = TRUE,
                                  scale = TRUE)

# Save matrix with all genes to file for downstream plotting
readr::write_rds(long_zscored_expression,
                 file.path(subset_dir, "hgg_zscored_expression.RDS"))

#### Filter focal CN data ------------------------------------------------------

# Filter focal CN to ATRT samples only
cn_metadata <- cn_df %>%
  dplyr::left_join(select_metadata,
                   by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::select(gene_symbol,
                sample_id,
                Kids_First_Participant_ID,
                biospecimen_id,
                status,
                cytoband) %>%
  dplyr::filter(sample_id %in% hgg_metadata_df$sample_id) %>%
  dplyr::distinct() # Remove duplicate rows produced as a result of not
                    # including the copy number variable from `cn_df`

# Write to file
readr::write_tsv(cn_metadata, file.path(subset_dir, "hgg_focal_cn.tsv.gz"))

#### Filter fusion data --------------------------------------------------------

fusion_df <- fusion_df %>%
  dplyr::select(Sample, FusionName) %>%
  dplyr::left_join(select_metadata,
                   by = c("Sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::filter(sample_id %in% hgg_metadata_df$sample_id)

# Write to file
readr::write_tsv(fusion_df, file.path(subset_dir, "hgg_fusion.tsv"))
