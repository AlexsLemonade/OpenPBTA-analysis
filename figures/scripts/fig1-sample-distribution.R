# JN Taroni and SJ Spielman for ALSF CCDL 2021-2022
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

#### Directories ---------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
figures_dir <- file.path(root_dir, "figures")
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")

# Tumor descriptor plots
main_output_dir <- file.path(figures_dir, "pdfs", "fig1", "panels")
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

# Define zenodo CSV output file
fig1b_csv <- file.path(zenodo_tables_dir, "figure-1b-data.csv")


#### Read in histologies -------------------------------------------------------

# Metadata we'll use to calculate counts
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"))
histologies_df <- histologies_df %>%
  filter(tumor_descriptor != "Unavailable")

#### Palette setup -------------------------------------------------------------

# Read in the cancer group palette
palette_df <- read_tsv(file.path(figures_dir,
                                 "palettes",
                                 "broad_histology_cancer_group_palette.tsv"))

# Tumor descriptor palette as a named vector
tumor_descriptor_palette <- read_tsv(file.path(figures_dir,
                                               "palettes",
                                               "tumor_descriptor_palette.tsv")) %>%
  tibble::deframe()

# Will mostly be using the cancer group display palette
cancer_group_palette <- palette_df$cancer_group_hex
names(cancer_group_palette) <- palette_df$cancer_group_display



#### Prep data for plotting ----------------------------------------------------

## Stacked bar plot that shows the assay types ##
experimental_strategy_df <- histologies_df %>%
  left_join(palette_df, by = c("broad_histology", "cancer_group")) %>%
  # Remove the Normal samples
  filter(sample_type != "Normal") %>%
  group_by(sample_id,  # We use sample_id to join genomic + transcriptomic samples
           broad_histology,
           cancer_group_display) %>%
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
           cancer_group_display,
           experimental_strategy) %>%
  count() %>%
  ungroup() %>%
  # Drop benign tumors and pre-cancerous lesions (i.e., NA in cancer_group)
  filter(!is.na(cancer_group_display)) %>%
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
  group_by(cancer_group_display) %>%
  mutate(cancer_group_n = sum(n))

# Set the y coordinate for the label based on the counts
cancer_group_counts_df <- experimental_strategy_df %>%
  select(broad_histology_display,
         cancer_group_display,
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
  select(-broad_histology) %>%
  ungroup()

# We'll order the panels by the display order in the palettes data frame
broad_histologies <- palette_df %>%
  select(broad_histology_display, broad_histology_order) %>%
  distinct() %>%
  arrange(broad_histology_order) %>%
  pull(broad_histology_display)

#### Tumor descriptor stacked bar plots ----------------------------------------

# a little more data prep -
data_descriptor_plot <- histologies_df %>%
  select(sample_id,
         broad_histology,
         cancer_group,
         tumor_descriptor) %>%
  # Drop benign tumors and pre-cancerous lesions (i.e., NA in cancer_group)
  drop_na(cancer_group) %>%
  # Join in display names and hex!
  inner_join(
    select(palette_df,
           broad_histology,
           broad_histology_hex,
           broad_histology_display,
           broad_histology_order),
    by = c("broad_histology")) %>%
  inner_join(
    select(palette_df,
           broad_histology,
           cancer_group,
           cancer_group_display,
           cancer_group_abbreviation),
    by = c("broad_histology", "cancer_group")) %>%
  distinct()

# The NAs in `cancer_group_abbreviation` should _all_ be associated with a cancer_display_group "Other" if we've joined up correctly -
abbr_check <- data_descriptor_plot %>%
  filter(is.na(cancer_group_abbreviation)) %>%
  pull(cancer_group_display)
if (!(all(abbr_check == "Other"))) stop("Wrangling bug setting cancer group abbreviations.")

# make the plot
data_descriptor_plot <- data_descriptor_plot %>%
  mutate(
    cancer_group_abbreviation = if_else(is.na(cancer_group_abbreviation), "Other", cancer_group_abbreviation),
    bh_strip = stringr::str_wrap(broad_histology_display, 18),
    bh_strip = forcats::fct_reorder(bh_strip, broad_histology_order),
    abbr = forcats::fct_relevel(cancer_group_abbreviation, "Other", after=Inf)
  )

descriptor_plot <- ggplot(data_descriptor_plot) +
  aes(x = abbr,
      fill = tumor_descriptor) +
  geom_bar(color = "black", size = 0.2) +
  scale_fill_manual(values = tumor_descriptor_palette) +
  facet_wrap(~bh_strip,
             scales = "free",
             nrow = 3) +
  labs(
    x = "Cancer group",
    y = "Number of tumors",
    fill = "Tumor descriptor"
  ) +
  ggpubr::theme_pubr() +
  guides(fill = guide_legend(nrow = 3)) +
  theme(axis.text.x = element_text(size = 7.5,
                                   angle = 45,
                                   hjust = 0.8,
                                   vjust = 0.9),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size = 0.3),
        # Single size since ggplot will not accept a per-panel vector here
        strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.4),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(0.01, 0.01, 0.0, 0.01, unit = "cm"),
        panel.spacing = unit(0.02, "cm"), 
        plot.margin = margin(0.01, 0.1, 0.01, 0.05, unit = "cm"))

#Save!
ggsave(filename = file.path(main_output_dir,
                            "tumor_descriptor_proportion_panel.pdf"),
       plot = descriptor_plot,
       width = 5.75,
       height = 6)

# Export figure data for zenodo upload
data_descriptor_plot %>%
  # remove columns with \n and other unneeded columns
  dplyr::select(-bh_strip, -abbr, -broad_histology, -broad_histology_hex, -broad_histology_order) %>%
  dplyr::arrange(sample_id) %>%
  readr::write_csv(fig1b_csv)

