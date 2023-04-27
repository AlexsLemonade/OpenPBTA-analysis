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

# Zenodo table directory
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")

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

# Zenodo data CSV
fig3d_csv <- file.path(zenodo_tables_dir, "figure-3d-data.csv")


#### Read in data --------------------------------------------------------------

# Read in clinical data
histologies_df <- read_tsv(histologies_file, guess_max = 10000)

# Read in chromothripsis per sample
chromo_per_sample_df <- read_tsv(chromo_file)

# Read in most recent cancer group palette
palette_df <- read_tsv(palette_file)

#### Prep ----------------------------------------------------------------------

# Need a data frame that contains the identifiers + cancer group
cancer_group_df <- histologies_df %>%
  filter(!is.na(pathology_diagnosis),
         !is.na(cancer_group)) %>%
  left_join(palette_df, by = c("broad_histology", "cancer_group")) %>%
  select(Kids_First_Biospecimen_ID, broad_histology, cancer_group_display, cancer_group, cancer_group_hex)

# Add cancer group and palette information
chromo_per_sample_df <- chromo_per_sample_df %>%
  left_join(cancer_group_df) %>%
  mutate(was_other = ifelse(cancer_group_display == "Other", "yes", "no"),
        cancer_group_display = case_when(cancer_group_display == "Other" ~ cancer_group,
                                        TRUE ~ cancer_group_display))


# We'll only be using the cancer group palette and we need a named vector of
# hexcodes
cancer_group_col_df <- chromo_per_sample_df %>%
  select(cancer_group_display, cancer_group_hex) %>%
  distinct()

cancer_group_palette <- cancer_group_col_df$cancer_group_hex
names(cancer_group_palette) <- cancer_group_col_df$cancer_group_display


### Check to determine that the correct Ns are being plotted with "Other" cancer_group_display --> cancer_group 
# Previously, we had 9/60 in the "Other" group:
previously_other <- chromo_per_sample_df %>%
  filter(was_other == "yes") %>%
  count(cancer_group_display)
sum(previously_other$n)

# N samples not being plotted because N < 3 is 22:
n_less_3 <- chromo_per_sample_df %>%
  filter(was_other == "yes") %>%
  count(cancer_group_display) %>%
  filter(n < 3)
sum(n_less_3$n)

# N samples being plotted because N >= 3 is 38:
n_greater_3 <- chromo_per_sample_df %>%
  filter(was_other == "yes") %>%
  count(cancer_group_display) %>%
  filter(n >=3)
sum(n_greater_3$n)


### Repeat all of these steps for cancer_group instead of display_group
chromoth_cancer_group_df <- chromo_per_sample_df %>%
  count(any_regions_logical, cancer_group_display) %>%
  tidyr::spread(key = any_regions_logical, value = n, fill=0) %>%
  group_by(cancer_group_display) %>%
  mutate(cancer_group_size = sum(`TRUE`, `FALSE`)) %>%
  mutate(prop = `TRUE` / cancer_group_size) %>%
  mutate(labels = paste0(`TRUE`, " / ", cancer_group_size)) %>%
  ungroup(cancer_group_display) %>%
  filter(cancer_group_size >= 3 & !is.na(cancer_group_display)) %>% # Only keep groups with >=3 tumors
  # Reorder cancer_group based on proportion
  mutate(cancer_group_display = fct_reorder(cancer_group_display, prop))

#### Plotting ------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_pubr())

# Save ggplot2 options
plot_options <- list(
  
)

bar_plot <- ggplot(chromoth_cancer_group_df ,
                   aes(x = cancer_group_display,
                       y = prop,
                       fill = cancer_group_display)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = labels), hjust = -0.2, size = 1.75) +
  scale_fill_manual(values = cancer_group_palette) +
  coord_flip() +
  # flipped, so thiss will apply to x:
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0, 0.01)) +
  labs(
    x = "",
    y = "Proportion of Tumors with\nChromothripsis Events"
  ) + 
  theme(
    axis.text.y = element_text(size = rel(0.5)), 
    axis.text.x = element_text(size = rel(0.475)), 
    axis.title = element_text(size = rel(0.5)),
    axis.line = element_line(size = rel(0.35)),
    axis.ticks = element_line(size = rel(0.35)), 
    plot.margin = margin(0,5,0,0)
  )

ggsave(file.path(output_dir, "chromothripsis_prop_cancer_group.pdf"),
       plot = bar_plot,
       height = 3,
       width = 3)

# Write CSV file for Zenodo upload, using the non-summarized data version
chromo_per_sample_df %>%
  dplyr::select(Kids_First_Biospecimen_ID, 
                # give this column an informative name
                chromothripsis_event_detected = any_regions_logical, 
                cancer_group_display) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  readr::write_csv(fig3d_csv)
