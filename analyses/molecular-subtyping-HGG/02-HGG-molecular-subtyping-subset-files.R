# This script subsets the focal copy number, RNA expression, fusion and
# histologies` and GISTIC's broad values files to include only High-grade glioma
# samples.

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
                Kids_First_Biospecimen_ID)

# Read in RNA expression data
stranded_expression <-
  readr::read_rds(
    file.path(
      root_dir,
      "data",
      "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
    )
  )

polya_expression <-
  readr::read_rds(
    file.path(
      root_dir,
      "data",
      "pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"
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

# Read in fusion data
fusion_df <- readr::read_tsv(
  file.path(root_dir, "data", "pbta-fusion-putative-oncogenic.tsv"))

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

# Read in snv consensus mutation data
snv_maf_df <-
  data.table::fread(file.path(root_dir,
                              "data",
                              "pbta-snv-consensus-mutation.maf.tsv.gz"),
                    select = c("Chromosome",
                               "Start_Position",
                               "End_Position",
                               "Strand",
                               "Variant_Classification",
                               "Tumor_Sample_Barcode",
                               "Hugo_Symbol",
                               "HGVSp_Short"),
                    data.table = FALSE)

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
# for samples classified as HGAT or those with defining lesions
hgg_lesions_df <- hgg_lesions_df %>%
  dplyr::filter(
    short_histology == "HGAT" |
      grepl("H3 G35 mutant|H3 K28 mutant", disease_type_reclassified)
  )

#### Filter metadata -----------------------------------------------------------

# Filter metadata for `High-grade glioma` and samples that should be classified
# as High-grade glioma based on defining lesions
hgg_metadata_df <- metadata %>%
  dplyr::filter(
    disease_type_new == "High-grade glioma" |
      sample_id %in% hgg_lesions_df$sample_id,
    sample_type == "Tumor",
    composition == "Solid Tissue"
  )

#### Filter expression data ----------------------------------------------------

filter_process_expression <- function(expression_mat) {
  # This function takes the collapsed FPKM expression matrix, selects relevant
  # columns (samples) via the Kids_First_Biospecimen_ID identifier, and then
  # log2(x + 1) transforms and z-scores the filtered matrix gene-wise.
  # It returns the z-scored matrix, where the columns are genes and the rows
  # are samples (biospecimen ID are the rownames).
  #
  # Only intended for use in the context of this script!

  # Filter to HGG samples only -- we can use hgg_metadata_df because it is
  # subset to RNA-seq samples
  filtered_expression <- expression_mat %>%
    dplyr::select(intersect(hgg_metadata_df$Kids_First_Biospecimen_ID,
                            colnames(expression_mat)))

  # Log2 transformation
  log_expression <- log2(filtered_expression + 1)

  # Scale does column centering, so we transpose first
  long_zscored_expression <- scale(t(log_expression),
                                   center = TRUE,
                                   scale = TRUE)
  return(long_zscored_expression)
}

# Save matrix with all genes to file for downstream plotting
filter_process_expression(stranded_expression) %>%
  readr::write_rds(file.path(subset_dir,
                             "hgg_zscored_expression.stranded.RDS"))
filter_process_expression(polya_expression) %>%
  readr::write_rds(file.path(subset_dir,
                             "hgg_zscored_expression.polya.RDS"))

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
  dplyr::filter(biospecimen_id %in% hgg_metadata_df$Kids_First_Biospecimen_ID) %>%
  dplyr::distinct() # Remove duplicate rows produced as a result of not
                    # including the copy number variable from `cn_df`

# Write to file
readr::write_tsv(cn_metadata, file.path(subset_dir, "hgg_focal_cn.tsv.gz"))

#### Filter fusion data --------------------------------------------------------

fusion_df <- fusion_df %>%
  dplyr::select(Sample, FusionName) %>%
  dplyr::left_join(select_metadata,
                   by = c("Sample" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::filter(Sample %in% hgg_metadata_df$Kids_First_Biospecimen_ID)

# Write to file
readr::write_tsv(fusion_df, file.path(subset_dir, "hgg_fusion.tsv"))

#### Filter GISTIC data --------------------------------------------------------

gistic_df <- gistic_df %>%
  tibble::column_to_rownames("Chromosome Arm") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::left_join(select_metadata, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% hgg_metadata_df$Kids_First_Biospecimen_ID) %>%
  dplyr::select(sample_id,
                Kids_First_Biospecimen_ID,
                `1p`,
                `19q`,
                `7p`,
                `7q`,
                `10p`,
                `10q`) # Select only the chromosome arms we are interested in

# Write to file
readr::write_tsv(gistic_df,
                 file.path(subset_dir, "hgg_gistic_broad_values.tsv"))

#### Filter SNV consensus maf data ---------------------------------------------

snv_maf_df <- snv_maf_df %>%
  dplyr::left_join(select_metadata,
                   by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::filter(Tumor_Sample_Barcode %in% hgg_metadata_df$Kids_First_Biospecimen_ID)

# Write to file
readr::write_tsv(snv_maf_df,
                 file.path(subset_dir, "hgg_snv_maf.tsv.gz"))
