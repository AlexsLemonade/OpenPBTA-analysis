# JN Taroni for ALSF CCDL 2021
#
# Create panels for representing sample distribution:
#   - Cancer group
#   - Experimental strategy
#   - Tumor distribution
#
# Each broad histology display group has an individual panel

#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggpattern)
library(patchwork)

#### Directories ---------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
figures_dir <- file.path(root_dir, "figures")
output_dir <- file.path(figures_dir, "pdfs", "fig1", "panels")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#### Read in histologies -------------------------------------------------------

# Metadata we'll use to calculate counts
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"))

#### Palette setup -------------------------------------------------------------

# Read in the cancer group palette
palette_df <- read_tsv(file.path(figures_dir,
                                 "palettes",
                                 "broad_histology_cancer_group_palette.tsv"))

# Will mostly be using the cancer group palette
cancer_group_palette <- palette_df$cancer_group_hex
names(cancer_group_palette) <- palette_df$cancer_group

# Tumor descriptor palette as a named vector
tumor_descriptor_palette <- read_tsv(file.path(figures_dir,
                                               "palettes",
                                               "tumor_descriptor_palette.tsv")) %>%
  tibble::deframe()

# "Palette" we'll use for the pattern in one of the stacked bar plots
experimental_strategy_pattern <- c(
  Both = "none",
  `RNA-Seq` = "crosshatch",
  `DNA-Seq` = "stripe"
)

#### Prep data for plotting ----------------------------------------------------

experimental_strategy_df <- histologies_df %>%
  # Remove the Normal samples
  filter(sample_type != "Normal") %>%
  group_by(sample_id,  # We use sample_id to join genomic + transcriptomic samples
           broad_histology,
           cancer_group) %>%
  # Summarize the experimental strategy -- we'll use this column to mutate in
  # the very next step
  summarize(summarized_experimental_strategy = paste(sort(unique(experimental_strategy)),
                                                     collapse = ", ")) %>%
  mutate(experimental_strategy = case_when(
    # When there's RNA-Seq and WGS|WXS|Targeted Sequencing
    str_detect(summarized_experimental_strategy, "RNA-Seq") &
      str_detect(summarized_experimental_strategy, "WGS|WXS|Targeted Sequencing") ~ "Both",
    # When there's RNA-Seq and no WGS|WXS|Targeted Sequencing
    str_detect(summarized_experimental_strategy, "RNA-Seq") &
      str_detect(summarized_experimental_strategy, "WGS|WXS|Targeted Sequencing", negate = TRUE) ~ "RNA-Seq",
    # When there's WGS|WXS|Targeted Sequencing and no RNA-seq
    str_detect(summarized_experimental_strategy, "RNA-Seq", negate = TRUE) &
      str_detect(summarized_experimental_strategy, "WGS|WXS|Targeted Sequencing") ~ "DNA-Seq"
  )) %>%
  # Drop the sample IDs
  ungroup() %>%
  select(-sample_id) %>%
  # New grouping to count
  group_by(broad_histology,
           cancer_group,
           experimental_strategy) %>%
  count() %>%
  ungroup() %>%
  # Drop benign tumors and pre-cancerous lesions (i.e., NA in cancer_group)
  filter(!is.na(cancer_group)) %>%
  # Reorder the experimental strategy for plotting
  mutate(experimental_strategy = factor(experimental_strategy,
                                        levels = c("Both",
                                                   "RNA-Seq",
                                                   "DNA-Seq"))) %>%
  inner_join(select(palette_df,
                    broad_histology,
                    broad_histology_display),
             by = "broad_histology") %>%
  distinct() %>%
  # Tack on the cancer group sample size
  group_by(cancer_group) %>%
  mutate(cancer_group_n = sum(n))

# Set the y coordinate for the label based on the counts
cancer_group_counts_df <- experimental_strategy_df %>%
  select(broad_histology_display,
         cancer_group,
         cancer_group_n) %>%
  distinct() %>%
  group_by(broad_histology_display) %>%
  mutate(y_coord = ceiling(cancer_group_n + (max(cancer_group_n) * 0.05)))

# We need to handle cases where there is only one cancer group in the broad
# histology display group otherwise the y coordinate for the label is way far
# out
broad_histology_counts <- cancer_group_counts_df %>%
  count(broad_histology_display, name = "broad_histology_n")

# And final data frame for plotting
experimental_strategy_df <- cancer_group_counts_df %>%
  # Add the broad histology counts required for the mutate step
  left_join(broad_histology_counts) %>%
  mutate(y_coord = case_when(
    # When there's only one group, we avoid doubling the y-coordinate
    broad_histology_n == 1 ~ (cancer_group_n * 1.05),
    TRUE ~ y_coord
  )) %>%
  select(-broad_histology_n) %>%
  left_join(experimental_strategy_df) %>%
  select(-broad_histology)

# Calculate the proportion of each tumor descriptor category for a cancer group
tumor_descriptor_df <- histologies_df %>%
  # Remove the Normal samples
  filter(sample_type != "Normal",
         !is.na(tumor_descriptor),
         !is.na(cancer_group)) %>%
  # Count needs to include the tumor_descriptor
  count(broad_histology,
        cancer_group,
        tumor_descriptor) %>%
  # Need to drop tumor_descriptor to calculate the proportion
  group_by(broad_histology,
           cancer_group) %>%
  # Calculate the proportin
  mutate(proportion = n / sum(n)) %>%
  inner_join(select(palette_df,
                    broad_histology,
                    broad_histology_display),
             by = "broad_histology") %>%
  distinct()

#### Custom function for plots by broad histology display ----------------------

# Function for creating the plots that show the sample size, experimental
# strategy, and tumor descriptor
create_two_plot_panel <- function(broad_histology_label)  {

  # Create a bar plot that shows the sample size, where the stacked pattern
  # shows the experimental strategy
  exp_strat_plot <- experimental_strategy_df %>%
    filter(broad_histology_display == broad_histology_label) %>%
    ggplot(aes(x = cancer_group,
               y = n)) +
    geom_bar_pattern(aes(fill = cancer_group,
                         pattern = experimental_strategy),
                     color = "#000000",
                     pattern_colour = "#000000",
                     pattern_fill = "#000000",
                     pattern_angle = 45,
                     pattern_spacing = 0.025,
                     stat = "identity") +
    geom_text(aes(y = y_coord,
                  label = cancer_group_n)) +
    scale_fill_manual(values = cancer_group_palette) +
    scale_pattern_manual(values = experimental_strategy_pattern) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 1,
                                     hjust = 1)) +
    labs(x = "",
         y = "Sample Size") +
    guides(pattern = FALSE,
           fill = FALSE,
           color = FALSE)

  # Create a stacked bar plot that shows the proportion of tumors in a cancer
  # group in each tumor descriptor category
  descriptor_plot <- tumor_descriptor_df %>%
    filter(broad_histology_display == broad_histology_label) %>%
    ggplot(aes(x = cancer_group,
               y = proportion,
               fill = tumor_descriptor)) +
    geom_bar(color = "#000000", position = "fill", stat = "identity") +
    scale_fill_manual(values = tumor_descriptor_palette) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 1,
                                     hjust = 1)) +
    labs(x = "",
         y = "Proportion") +
    guides(fill = FALSE)

  # Use patchwork to put plots side by side
  two_plot_panel <- exp_strat_plot + descriptor_plot +
    plot_annotation(title = broad_histology_label) &
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  if (broad_histology_label == "Other") {
    # Add a unified x-axis label because this would be the bottom plot
    two_plot_panel <- gridExtra::grid.arrange(patchworkGrob(two_plot_panel),
                                              bottom = "Cancer Group")
  }

  # Save PDF of two plot panel
  ggsave(file.path(output_dir,
                   paste(broad_histology_label, "two_panels.pdf")),
         plot = two_plot_panel,
         width = 7,
         height = 7)

}

#### Plots by broad histology --------------------------------------------------

# We'll order the panels by the display order in the palettes data frame
broad_histologies <- palette_df %>%
  select(broad_histology_display, broad_histology_order) %>%
  distinct() %>%
  arrange(broad_histology_order) %>%
  pull(broad_histology_display)

# Create PDFs of two plot panels - one per broad histology display group
purrr::walk(broad_histologies,
            ~ create_two_plot_panel(.x))

#### Legends -------------------------------------------------------------------
# We pick Diffuse astrocytic and oligodendroglial tumor because every field is
# represented

# Save the legend for the patterns for experimental strategy
as_ggplot(get_legend(
  experimental_strategy_df %>%
    filter(broad_histology_display == "Diffuse astrocytic and oligodendroglial tumor") %>%
    ggplot(aes(x = cancer_group,
               y = n)) +
    geom_bar_pattern(aes(fill = cancer_group,
                         pattern = experimental_strategy),
                     color = "#000000",
                     pattern_colour = "#000000",
                     pattern_fill = "#000000",
                     pattern_angle = 45,
                     pattern_spacing = 0.025,
                     stat = "identity") +
    geom_text(aes(y = y_coord,
                  label = cancer_group_n)) +
    scale_fill_manual(values = cancer_group_palette) +
    scale_pattern_manual(values = experimental_strategy_pattern) +
    theme_classic() +
    labs(pattern = "Assay") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = FALSE)
)) %>%
  ggsave(filename = file.path(output_dir, "assay_pattern_legend.pdf"),
         height = 2,
         width = 2)

# Save the legend for the tumor descriptors
as_ggplot(get_legend(
  tumor_descriptor_df %>%
    filter(broad_histology_display == "Diffuse astrocytic and oligodendroglial tumor") %>%
    ggplot(aes(x = cancer_group,
               y = proportion,
               fill = tumor_descriptor)) +
    geom_bar(color = "#000000", position = "fill", stat = "identity") +
    scale_fill_manual(values = tumor_descriptor_palette) +
    theme_classic() +
    labs(fill = "Tumor Type")
)) %>%
  ggsave(filename = file.path(output_dir, "tumor_descriptor_legend.pdf"),
         height = 3,
         width = 4)
