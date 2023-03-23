# S. Spielman for ALSF CCDL, 2022
#
# This script makes panel 4E, the hypermutators heatmap, and its legend
# Adapted from `mutational-signatures/08-explore_hypermutators.Rmd


library(ComplexHeatmap)
library(tidyverse)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig4", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig4e_csv <- file.path(zenodo_tables_dir, "figure-4e-data.csv")

# Directory with result data to input for plot
input_dir <- file.path(root_dir, 
                       "analyses", 
                       "mutational-signatures", 
                       "results")

## Input and output files -------------------------------------
meta_file    <- file.path(root_dir, "data", "pbta-histologies.tsv")
palette_file <- file.path(root_dir, "figures", "palettes", "broad_histology_cancer_group_palette.tsv")
bin_palette_file <- file.path(root_dir, "figures", "palettes", "binary_color_palette.tsv")
td_palette_file <- file.path(root_dir, "figures", "palettes", "tumor_descriptor_palette.tsv")

tmb_coding_file <- file.path(root_dir, "data", "pbta-snv-consensus-mutation-tmb-coding.tsv")
signatures_file  <- file.path(input_dir, "deconstructsigs_exposures_merged.tsv")

heatmap_pdf <- file.path(output_dir, "hypermutator_sigs_heatmap.pdf")
heatmap_legends_pdf <- file.path(output_dir, "hypermutator_sigs_heatmap_legends.pdf")


## Prepare data for plotting --------------------------------------------------

# Read in tumor descriptor and palettes
td_palette <- readr::read_tsv(td_palette_file)
bin_palette <- readr::read_tsv(bin_palette_file) %>%
  mutate(color_names = case_when(
    color_names == "binary_1" ~ "Solid Tissue",
    color_names == "binary_2" ~ "Derived Cell Line")) %>%
  filter(!is.na(color_names))

# Read in metadata
meta <- readr::read_tsv(meta_file, guess_max = 10000)

# Read in exposures
exposures <- readr::read_tsv(signatures_file)

# Read in coding tmb data
tmb_coding <- readr::read_tsv(tmb_coding_file)



## Make heatmap ------------------------

# Find all hypermutators
hypermutator_df <- tmb_coding %>%
  # filter to the hypermutators
  filter(tmb >= 10) %>%
  # keep only tmb and ID columns
  select(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode, tmb)
 
# Prepare signatures to merge in
sigs <- exposures %>%
  filter(signature != "Other") %>%
  spread(signature, exposure) %>%
  select(Kids_First_Biospecimen_ID, `1`:N6)

# Create signature matrix for only patients with hypermutant samples
# Gather all BS_IDs (soi = "samples of interest"), whether hypermutant or not, for pts with hypermutant samples
soi_df <- meta %>%
  filter(Kids_First_Biospecimen_ID %in% hypermutator_df$Kids_First_Biospecimen_ID) %>%
  select(Kids_First_Participant_ID) %>%
  left_join(meta) %>%
  #keep only WGS samples; remove normals; keep cell lines
  filter(experimental_strategy == "WGS" & !is.na(pathology_diagnosis)) %>%
  select(Kids_First_Participant_ID, 
         Kids_First_Biospecimen_ID, 
         tumor_descriptor, 
         composition) %>%
  distinct() %>%
  inner_join(sigs) %>%
  # Arrange here on participant frequency to set an overall heatmap order
  mutate(Kids_First_Participant_ID = fct_infreq(Kids_First_Participant_ID)) %>%
  arrange(Kids_First_Participant_ID)

# Use this order moving forward to ensure legend matches heatmap patient order
soi_order <- soi_df$Kids_First_Biospecimen_ID

# Filter sigs to only relevant samples
sigs_soi <- sigs %>% 
  filter(Kids_First_Biospecimen_ID %in% soi_df$Kids_First_Biospecimen_ID) %>%
  mutate(Kids_First_Biospecimen_ID = fct_relevel(Kids_First_Biospecimen_ID, soi_order)) %>%
  arrange(Kids_First_Biospecimen_ID) %>%
  as.data.frame()

# add spaces at the end of BS_ids to make a space between the heatmap and legend
rownames(sigs_soi) <- paste0(soi_order, "     ")

sigs_soi$Kids_First_Biospecimen_ID <- NULL
sigs_soi_matrix <- as.matrix(sigs_soi)

# Create annotation dataframe
anno <- tmb_coding %>%
  mutate(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  select(Kids_First_Biospecimen_ID, tmb) %>%
  # add a column for mutation status
  mutate(
    `Mutation status` = case_when(
      tmb < 10 ~ "Normal",
      tmb >= 10 & tmb < 100 ~ "Hypermutant",
      tmb >= 100 ~ "Ultra-hypermutant")
  ) %>%
  right_join(soi_df, by = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  # mutate names for annotation
  mutate(Patient = Kids_First_Participant_ID,
         `Phase of Therapy` = tumor_descriptor,
         Composition = composition) %>%
  select(Patient, Kids_First_Biospecimen_ID, `Phase of Therapy`, Composition, `Mutation status`) %>%
  mutate(Kids_First_Biospecimen_ID = fct_relevel(Kids_First_Biospecimen_ID, soi_order)) %>%
  arrange(Kids_First_Biospecimen_ID) %>%
  as.data.frame() 

# set rownames to match matrix, remove BS_ID
rownames(anno) <- paste0(anno$Kids_First_Biospecimen_ID, "     ")
anno$Kids_First_Biospecimen_ID <- NULL

# Add tumor descriptor palette from dataframe
td_palette_filtered <- td_palette %>%
  filter(color_names %in% anno$`Phase of Therapy`)
td_palette_col <- as.character(td_palette_filtered$hex_codes)
names(td_palette_col) <- td_palette_filtered$color_names

pt_palette <- colorRampPalette(colorblindr::palette_OkabeIto)(length(unique(anno$Patient)))
names(pt_palette) <- unique(anno$Patient)

comp_palette <- as.character(bin_palette$hex_codes)
names(comp_palette) <- unique(anno$Composition)

# Specify colors
ann_colors = list(
  # mutator color from [this PR](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/1280/files#diff-bc8134008706f38e7795b59696928058d123e44b4172df68105a26ad8aa0a398R298-R301) 
  `Mutation status` = c(Normal = "grey40",
                        Hypermutant = "orange",
                        `Ultra-hypermutant` = "red"),
  Composition = comp_palette,
  `Phase of Therapy` = td_palette_col,
  Patient = pt_palette)

# Heatmap annotation
row_anno = rowAnnotation(df = anno,
                         col = ann_colors, show_legend=FALSE)

# Make heatmap without legends
heat_plot <- Heatmap(sigs_soi_matrix[, c("1", "3", "8", "11", "18", "19", "MMR2", "N6")], 
                     name = "Signature weights",
                     col = colorRampPalette(c("#f1f1f1", "#2166ac"))(50),
                     cluster_rows = FALSE,
                     show_row_names = TRUE,
                     show_heatmap_legend=FALSE,
                     cluster_columns = FALSE,
                     right_annotation = row_anno,
                     rect_gp = gpar(col = "white"),
                     row_split = anno$Patient,
                     row_title = NULL, 
                     column_title = "RefSig Mutational Signature Weights", 
                     column_title_side = "bottom")

heat_plot

# Make separate legends to compile in Illustrator
patient_legend <- Legend(
  labels = names(ann_colors$Patient),
  legend_gp = gpar(fill = ann_colors$Patient),
  title = "Patient"
)

mutator_legend <- Legend(
  labels = names(ann_colors$`Mutation status`),
  legend_gp = gpar(fill = ann_colors$`Mutation status`),
  title = "Mutation status"
)

phase_legend <- Legend(
  labels = names(ann_colors$`Phase of Therapy`),
  legend_gp = gpar(fill = ann_colors$`Phase of Therapy`),
  title = "Phase of Therapy"
)

comp_legend <- Legend(
  labels = names(ann_colors$Composition),
  legend_gp = gpar(fill = ann_colors$Composition),
  title = "Composition"
)

weights_legend <- color_mapping_legend(heat_plot@matrix_color_mapping, plot = FALSE, 
                                       legend_direction = "horizontal")

heat_legends <- packLegend(patient_legend,
                           mutator_legend,
                           phase_legend,
                           comp_legend,
                           weights_legend,
                           direction = "horizontal",
                           column_gap = unit(0.75, "cm")
) 


# save heatmap 
pdf(heatmap_pdf, width = 8, height = 4)
print(heat_plot)
dev.off()

# save heatmap legends
pdf(heatmap_legends_pdf, width = 8, height = 2)
draw(heat_legends)
dev.off()  


# Export Figure 4E data in CSV format for Zenodo upload
soi_df %>%
  # arrange on partipant, which is shown in the MS figure
  dplyr::arrange(Kids_First_Participant_ID) %>%
  # export
  readr::write_csv(fig4e_csv)


