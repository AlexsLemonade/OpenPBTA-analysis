# JN Taroni for ALSF CCDL 2021
# Adapted from Laura Egolf (analyses/chromothripsis/03-plot-chromothripsis-by-histology.Rmd)
#
# Create a panel with a barplot counting the number of samples with
# chromothripsis per cancer group

library(tidyverse)
library(ggpubr)

#### Directories ---------------------------------------------------------------

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Figures directory
figures_dir <- file.path(root_dir, "figures")

# Declare output directory
output_dir <- file.path(figures_dir, "pdfs", "fig3", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directory
analyses_dir <- file.path(root_dir, "analyses")

# Chromothripsis results directory
chromothripsis_dir <- file.path(analyses_dir, "chromothripsis", "results")

#### Files ---------------------------------------------------------------------

# Histologies file
histologies_file <- file.path(data_dir, "pbta-histologies.tsv")

# Results file with chromothripsis info for each sample
chromo_file <- file.path(chromothripsis_dir,
                         "chromothripsis_summary_per_sample.txt")

# Most recent palette for cancer group
palette_file <- file.path(figures_dir,
                          "palettes",
                          "broad_histology_cancer_group_palette.tsv")

#### Read in data --------------------------------------------------------------

# Read in clinical data
histologies_df <- read_tsv(histologies_file, guess_max = 10000)

# Read in chromothripsis per sample
chromo_per_sample_df <- read_tsv(chromo_file)

# Read in most recent cancer group palette
palette_df <- read_tsv(palette_file)

#### Prep ----------------------------------------------------------------------

# We'll only be using the cancer group palette and we need a named vector of
# hexcodes
cancer_group_col_df <- palette_df %>%
  select(cancer_group, cancer_group_hex) %>%
  distinct()

# Can't do pull with names with our versions
cancer_group_palette <- cancer_group_col_df$cancer_group_hex
names(cancer_group_palette) <- cancer_group_col_df$cancer_group

# Need a data frame that contains the identifiers + cancer group
cancer_group_df <- histologies_df %>%
  select(Kids_First_Biospecimen_ID, cancer_group)

# Add cancer group and palette information t
chromo_per_sample_df <- chromo_per_sample_df %>%
  left_join(cancer_group_df)

### Repeat all of these steps for cancer_group instead of display_group
chromoth_cancer_group_df <- chromo_per_sample_df %>%
  count(any_regions_logical, cancer_group) %>%
  tidyr::spread(key = any_regions_logical, value = n, fill=0) %>%
  group_by(cancer_group) %>%
  mutate(cancer_group_size = sum(`TRUE`, `FALSE`)) %>%
  mutate(prop = `TRUE` / cancer_group_size) %>%
  mutate(labels = paste0(`TRUE`, " / ", cancer_group_size)) %>%
  ungroup(cancer_group) %>%
  filter(cancer_group_size >= 3 & !is.na(cancer_group)) %>% # Only keep groups with >=3 tumors
  # Reorder cancer_group based on proportion
  mutate(cancer_group = fct_reorder(cancer_group, prop))

#### Plotting ------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_pubr())

# Save ggplot2 options
plot_options <- list(
  ylim(c(0, 1)),
  xlab(NULL),
  ylab("Proportion of Tumors with Chromothripsis Events"),
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95))
)

bar_plot <- ggplot(chromoth_cancer_group_df ,
                   aes(x = cancer_group,
                       y = prop,
                       fill = cancer_group)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = labels), hjust = -0.2, size = 3.25) +
  scale_fill_manual(values = cancer_group_palette) +
  coord_flip() +
  plot_options

ggsave(file.path(output_dir, "chromothripsis_prop_cancer_group.pdf"),
       plot = bar_plot,
       height = 7,
       width = 9)
