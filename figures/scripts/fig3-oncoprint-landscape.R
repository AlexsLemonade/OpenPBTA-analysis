# Oncoprint Landscape Figure
#
# 2020
# Chante Bethell for ALSF - CCDL
#
# This script is intended to run steps needed to create Figure 3.

#### Set Up --------------------------------------------------------------------

# Load maftools
library(maftools)

# Load in patchwork for assembling the final multipanel figure
library(patchwork)

# Get `magrittr` pipe
`%>%` <- dplyr::`%>%`

# Define clinical features object for plotting -- this will be the same
# across maf objects
clinical_features = c(
  "broad_histology",
  "short_histology",
  "reported_gender",
  "tumor_descriptor",
  "molecular_subtype"
)

#### Directories and Files -----------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Path to the data obtained via `00-map-to-sample_id.R`
scratch_oncoprint_dir <- file.path(root_dir, "scratch", "oncoprint_files")

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

# Path to oncoprint landscape directory
oncoprint_dir <- file.path(root_dir, "analyses", "oncoprint-landscape")

# Source the color palette for the oncoprint plots
source(
  file.path(
    oncoprint_dir,
    "util",
    "oncoplot-palette.R"
  )
)

#### Read in data --------------------------------------------------------------

# Read in metadata
metadata <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv")) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_id)

# Read in MAF files
maf_df_primary_plus <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary-plus_maf.tsv"
  ))

maf_df_primary_only <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary_only_maf.tsv"
  ))

# Read in cnv files
cnv_df_primary_plus <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary-plus_cnv.tsv"
  ))

cnv_df_primary_only <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary_only_cnv.tsv"
  ))

# Read in fusion files
fusion_df_primary_plus <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary-plus_fusions.tsv"
  ))

fusion_df_primary_only <-
  readr::read_tsv(file.path(
    scratch_oncoprint_dir,
    "all_participants_primary_only_fusions.tsv"
  ))

# Read in recurrent focal CNVs file
recurrent_focal_cnvs <-
  readr::read_tsv(
    file.path(
      root_dir,
      "analyses",
      "focal-cn-file-preparation",
      "results",
      "consensus_seg_recurrent_focal_cn_units.tsv"
    )
  )

#### Functions ----------------------------------------------------------

source(file.path(
  oncoprint_dir,
  "util",
  "oncoplot-functions.R"
))


#### Format recurrent focal CN objects ----------------------------------------

# Let's filter the recurrent focal calls data frame to include the primary
# plus samples
most_focal_recurrent_primary_plus <- recurrent_focal_cnvs %>%
  filter(status != "uncallable") %>%
  # Join the metadata to get the `Tumor_Sample_Barcode` column
  left_join(metadata, by = "Kids_First_Biospecimen_ID") %>%
  # Select and rename the needed columns for creating the maf object
  select(Hugo_Symbol = region, Tumor_Sample_Barcode, Variant_Classification = status) %>%
  # Filter for primary plus samples
  filter(Tumor_Sample_Barcode %in% cnv_df_primary_plus$Tumor_Sample_Barcode)

# Filter the recurrent focal calls data frame to include the primary only
# samples
most_focal_recurrent_primary_only <- recurrent_focal_cnvs %>%
  filter(status != "uncallable") %>%
  # Join the metadata to get the `Tumor_Sample_Barcode` column
  left_join(metadata, by = "Kids_First_Biospecimen_ID") %>%
  # Select and rename the needed columns for creating the maf object
  select(Hugo_Symbol = region, Tumor_Sample_Barcode, Variant_Classification = status) %>%
  # Filter for primary only samples
  filter(Tumor_Sample_Barcode %in% cnv_df_primary_only$Tumor_Sample_Barcode)

#### Generate MAF objects for plotting ----------------------------------------

# Generate the primary-plus maf_object
primary_plus <- prepare_maf_object(
  maf_df = maf_df_primary_plus,
  cnv_df = most_focal_recurrent_primary_plus,
  metadata = metadata,
  fusion_df = fusion_df_primary_plus
)

# Generate the primary-only maf object
primary_only <- prepare_maf_object(
  maf_df = maf_df_primary_only,
  cnv_df = most_focal_recurrent_primary_only,
  metadata = metadata,
  fusion_df = fusion_df_primary_only
)

#### Generate oncoprint plots -------------------------------------------------

# Given a the `primary_plus` maf object, plot an oncoprint of the variants in
# the dataset and save as a png file in the scratch directory.
png(
  file.path(scratch_oncoprint_dir, "fig3a_oncoprint.png"),
  width = 45,
  height = 35,
  units = "cm",
  res = 300
)
oncoplot(
  primary_plus,
  clinicalFeatures = clinical_features,
  logColBar = TRUE,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = TRUE,
  removeNonMutated = TRUE,
  annotationFontSize = 0.8,
  SampleNamefontSize = 0.5,
  fontSize = 0.7,
  colors = color_palette,
  showTitle = FALSE
)
dev.off()

# Given a the `primary_only` maf object, plot an oncoprint of the variants in
# the dataset and save as a png file in the scratch directory.
png(
  file.path(scratch_oncoprint_dir, "fig3b_oncoprint.png"),
  width = 45,
  height = 35,
  units = "cm",
  res = 300
)
oncoplot(
  primary_only,
  clinicalFeatures = clinical_features,
  logColBar = TRUE,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = TRUE,
  removeNonMutated = TRUE,
  annotationFontSize = 0.8,
  SampleNamefontSize = 0.5,
  fontSize = 0.7,
  colors = color_palette,
  showTitle = FALSE
)
dev.off()

#### Assemble multipanel plot -------------------------------------------------

# Read in the `primary_plus` and `primary_only` oncoprints
primary_plus_onco <- magick::image_read(file.path(scratch_oncoprint_dir, "fig3a_oncoprint.png")) %>%
  magick::image_annotate( "Primary and Secondary Independent Samples Oncoprint", size = 70, gravity = "north", color = "black")
primary_only_onco <- magick::image_read(file.path(scratch_oncoprint_dir, "fig3b_oncoprint.png")) %>%
  magick::image_annotate("Primary Only Independent Samples Oncoprint", size = 70, gravity = "north", color = "black")

# Combine the oncoprints into one plot
combined_oncoprint <- magick::image_append(c(primary_plus_onco, primary_only_onco))

# Write the final image as a PNG file
magick::image_write(combined_oncoprint, file.path(output_dir, "fig3-oncoprint-lanscape.png"))
