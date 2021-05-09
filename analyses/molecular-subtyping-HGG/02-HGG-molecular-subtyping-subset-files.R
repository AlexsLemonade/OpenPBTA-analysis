# This script subsets the focal copy number, RNA expression, fusion and
# histologies` and GISTIC's broad values files to include only High-grade glioma
# samples.
#
# Chante Bethell for CCDL 2020
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript 'analyses/molecular-subtyping-ATRT/02-HGG-molecular-subtyping-subset-files.R'


#### Set Up --------------------------------------------------------------------

library(tidyverse)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
scratch_dir <- file.path(root_dir, "scratch")

# Set path to subset directory
subset_dir <-
  file.path(root_dir, "analyses", "molecular-subtyping-HGG", "hgg-subset")

if (!dir.exists(subset_dir)) {
  dir.create(subset_dir)
}

# Read in metadata
metadata <-
  read_tsv(file.path(root_dir, "data", "pbta-histologies-base.tsv"),
           guess_max = 10000)

# Select wanted columns in metadata for merging and assign to a new object
select_metadata <- metadata %>%
  select(sample_id,
         Kids_First_Participant_ID,
         Kids_First_Biospecimen_ID)

# Read in RNA expression data
stranded_expression <-
  read_rds(
    file.path(
      root_dir,
      "analyses",
      "collapse-rnaseq",
      "results",
      "pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
    )
  )

polya_expression <-
  read_rds(
    file.path(
      root_dir,
      "analyses",
      "collapse-rnaseq",
      "results",
      "pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"
    )
  )


# Read in focal CN data
## TODO: If annotated files get included in data download
cn_df <- read_tsv(file.path(
  root_dir,
  "data",
  "consensus_seg_annotated_cn_autosomes.tsv.gz"
))

# Read in fusion data
fusion_df <- read_tsv(
  file.path(root_dir, "analyses","fusion_filtering", "results", "pbta-fusion-putative-oncogenic.tsv"))

# Read in GISTIC `broad_values_by_arm.txt` file
unzip(file.path(root_dir, "data", "pbta-cnv-consensus-gistic.zip"),
      exdir = file.path(root_dir, "data"),
      files = file.path("pbta-cnv-consensus-gistic", "broad_values_by_arm.txt"))

gistic_df <- data.table::fread(file.path(root_dir,
                                         "data",
                                         "pbta-cnv-consensus-gistic",
                                         "broad_values_by_arm.txt"),
                               data.table = FALSE)


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
hgg_lesions_df <- read_tsv(
  file.path(
    root_dir,
    "analyses",
    "molecular-subtyping-HGG",
    "results",
    "HGG_defining_lesions.tsv"
  )
)

# Read in the JSON file that contains the strings we'll use to include or
# exclude samples for subtyping - see 00-HGG-select-pathology-dx
path_dx_list <- jsonlite::fromJSON(
  file.path(subset_dir,
            "hgg_subtyping_path_dx_strings.json")
)

#### Filter metadata -----------------------------------------------------------

# Filter metadata based on pathology diagnosis fields and include samples that
# should be classified as high-grade glioma based on defining lesions

# First, filter to exclude cell lines
tumor_metadata_df <- metadata %>%
  filter(
    sample_type == "Tumor"
  )

# Samples included on the basis of the pathology diagnosis fields
path_dx_df <- tumor_metadata_df %>%
  # Inclusion on the basis of CBTTC harmonized pathology diagnoses
  filter(pathology_diagnosis %in% path_dx_list$exact_path_dx |
         # Inclusion based on pathology free text diagnosis
         pathology_free_text_diagnosis ==path_dx_list$gliomatosis_path_free_text_exact)
  

# Now samples on the basis of the defining lesions
hgg_sample_ids <- hgg_lesions_df %>%
  filter(defining_lesion) %>%
  pull(sample_id)
lesions_df <- tumor_metadata_df %>%
  filter(sample_id %in% hgg_sample_ids)

# Putting it all together now
hgg_metadata_df <- bind_rows(
  path_dx_df,
  lesions_df
) %>%
  # Remove duplicates
  distinct()

# Add a TSV that's the metadata for the samples that will be included in
# the subtyping
# Useful intermediate file for examination, but also can be used in the next
# notebook for inclusion criteria without repeating the logic
write_tsv(hgg_metadata_df,
          file.path(subset_dir, "hgg_metadata.tsv"))

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
    select(intersect(hgg_metadata_df$Kids_First_Biospecimen_ID,
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
  write_rds(file.path(subset_dir,
                      "hgg_zscored_expression.stranded.RDS"))
filter_process_expression(polya_expression) %>%
  write_rds(file.path(subset_dir,
                      "hgg_zscored_expression.polya.RDS"))

#### Filter focal CN data ------------------------------------------------------

# Filter focal CN to ATRT samples only
cn_metadata <- cn_df %>%
  left_join(select_metadata,
            by = c("biospecimen_id" = "Kids_First_Biospecimen_ID")) %>%
  select(gene_symbol,
         sample_id,
         Kids_First_Participant_ID,
         biospecimen_id,
         status,
         cytoband) %>%
  filter(biospecimen_id %in% hgg_metadata_df$Kids_First_Biospecimen_ID) %>%
  distinct() %>% # Remove duplicate rows produced as a result of not
  # including the copy number variable from `cn_df`
  arrange(Kids_First_Participant_ID, sample_id)

# Write to file
write_tsv(cn_metadata, file.path(subset_dir, "hgg_focal_cn.tsv.gz"))

#### Filter fusion data --------------------------------------------------------

fusion_df <- fusion_df %>%
  select(Sample, FusionName) %>%
  left_join(select_metadata,
            by = c("Sample" = "Kids_First_Biospecimen_ID")) %>%
  filter(Sample %in% hgg_metadata_df$Kids_First_Biospecimen_ID) %>%
  arrange(Sample, FusionName)

# Write to file
write_tsv(fusion_df, file.path(subset_dir, "hgg_fusion.tsv"))

#### Filter GISTIC data --------------------------------------------------------

gistic_df <- gistic_df %>%
  tibble::column_to_rownames("Chromosome Arm") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  left_join(select_metadata, by = "Kids_First_Biospecimen_ID") %>%
  filter(Kids_First_Biospecimen_ID %in% hgg_metadata_df$Kids_First_Biospecimen_ID) %>%
  arrange(Kids_First_Participant_ID, sample_id) %>%
  select(sample_id,
         Kids_First_Biospecimen_ID,
         `1p`,
         `19q`,
         `7p`,
         `7q`,
         `10p`,
         `10q`) # Select only the chromosome arms we are interested in

# Write to file
write_tsv(gistic_df,
          file.path(subset_dir, "hgg_gistic_broad_values.tsv"))

#### Filter SNV consensus maf data ---------------------------------------------

snv_maf_df <- snv_maf_df %>%
  left_join(select_metadata,
            by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  filter(Tumor_Sample_Barcode %in% hgg_metadata_df$Kids_First_Biospecimen_ID) %>%
  arrange(Kids_First_Participant_ID, sample_id)

# Write to file
write_tsv(snv_maf_df,
          file.path(subset_dir, "hgg_snv_maf.tsv.gz"))
