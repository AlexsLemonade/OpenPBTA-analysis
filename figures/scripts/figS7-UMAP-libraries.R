# S. Spielman for ALSF CCDL 2023
#
# This script creates a panel for Figure S7 of UMAP for
#  samples derived from _both_ polyA and stranded RNA-Seq library strategies.
# This script follows the approach in:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/selection-strategy-comparison/01-selection-strategies.rmd

library(tidyverse)
set.seed(2023) # umap seed


# Directories and files ---------------------------

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "supp", "figs7", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory and data files: metadata and RSEM files
data_dir <- file.path(root_dir, "data")
metadata_file <- file.path(data_dir, "pbta-histologies.tsv")
polya_file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm.polya.rds")
stranded_file <- file.path(data_dir, "pbta-gene-expression-rsem-fpkm.stranded.rds")

# Palette file
palette_file <- file.path(root_dir,
                          "figures",
                          "palettes",
                          "broad_histology_cancer_group_palette.tsv")

# Output PDF filename
output_pdf <- file.path(output_dir, "umap_polya_stranded.pdf")


# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
figS7b_csv <- file.path(zenodo_tables_dir, "figure-S7b-data.csv")

# Read and pre-process data -------------------------

# Read in metadata associated palette
metadata_df <- read_tsv(metadata_file, guess_max = 10000)
palette_df <- read_tsv(palette_file)


# Read in FPKM files, join files, and convert to matrix ------------------
# Derived from:
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/selection-strategy-comparison/01-selection-strategies.rmd#L127

# read in the data frames
polya <- read_rds(polya_file)
stranded <- read_rds(stranded_file)

# join polyA and stranded data together
exp_rsem <- bind_cols(polya, stranded[,-1]) %>%
  filter(complete.cases(.))

# transpose the expression values to a matrix
exp_rsem_t <- t(exp_rsem[,-1])
colnames(exp_rsem_t) <- exp_rsem[[1]]


# Perform UMAP ---------------------

# Remove low counts:
dm_set <- exp_rsem_t[, colSums(exp_rsem_t) > 100]


# Using 15 nn: https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/627ec427ad0a8d9d913e614c9db50546c56d8283/analyses/selection-strategy-comparison/01-selection-strategies.rmd#L10
combined_umap <- umap::umap(dm_set, n_neighbors = 15)

# Plot UMAP ------------------------

umap_df <- data.frame(combined_umap$layout) %>%
  rename_all(.funs = gsub, pattern = "^X", replacement = "UMAP") %>%
  rownames_to_column(var = "Kids_First_Biospecimen_ID") %>%
  # join with library strategy
  left_join(
    select(metadata_df,
           Kids_First_Biospecimen_ID,
           RNA_library,
           cancer_group),
    by = "Kids_First_Biospecimen_ID"
  ) %>%
  # join with cancer_group_display
  inner_join(
    select(palette_df,
           cancer_group,
           cancer_group_display),
    by = "cancer_group"
  ) %>%
  # remove NA cancer group displays
  drop_na(cancer_group_display)

# Set up cancer group display palette
umap_palette <- palette_df$cancer_group_hex
names(umap_palette) <- palette_df$cancer_group_display

umap_plot <- ggplot(umap_df) +
  aes(x = UMAP1,
      y = UMAP2,
      shape = RNA_library,
      color = cancer_group_display) +
  geom_point(alpha = 0.35,
             size = 1.5) +
  scale_shape_manual(name = "RNA library strategy", values = c(17, 19)) +
  scale_color_manual(name = "Cancer group", values = umap_palette) +
  ggpubr::theme_pubr() +
  theme(
    legend.direction = "vertical",
    legend.title = element_text(size = rel(0.6)),
    legend.text = element_text(size = rel(0.55)),
    legend.key.size = unit(0.3, "cm")
  )

ggsave(
  output_pdf,
  umap_plot,
  width = 7,
  height = 6
)

# Export CSV for Zenodo upload
umap_df %>%
  dplyr::select(Kids_First_Biospecimen_ID, everything()) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(figS7b_csv)




