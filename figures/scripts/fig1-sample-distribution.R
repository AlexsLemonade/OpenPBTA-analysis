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
library(ggpattern)

#### Directories ---------------------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
figures_dir <- file.path(root_dir, "figures")

# Tumor descriptor plots
main_output_dir <- file.path(figures_dir, "pdfs", "fig1", "panels")
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

# Assay type stacked bar plots
supp_output_dir <- file.path(figures_dir, "pdfs", "supp", "sample-distribution")
dir.create(supp_output_dir, recursive = TRUE, showWarnings = FALSE)

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
    by = "broad_histology") %>%
  inner_join(
    select(palette_df,
           cancer_group,
           cancer_group_display,
           cancer_group_abbreviation),
    by = "cancer_group") %>%
  distinct() 

# The NAs in `cancer_group_abbreviation` should _all_ be associated with a cancer_display_group "Other" if we've joined up correctly - 
abbr_check <- data_descriptor_plot %>%
  filter(is.na(cancer_group_abbreviation)) %>%
  pull(cancer_group_display)
if (!(all(abbr_check == "Other"))) stop("Wrangling bug setting cancer group abbreviations.")

#make the plot
descriptor_plot <- data_descriptor_plot %>%
  mutate(
    cancer_group_abbreviation = if_else(is.na(cancer_group_abbreviation), "Other", cancer_group_abbreviation),
    bh_strip = stringr::str_wrap(broad_histology_display, 25),
    bh_strip = forcats::fct_reorder(bh_strip, broad_histology_order),
    abbr = forcats::fct_relevel(cancer_group_abbreviation, "Other", after=Inf)
  ) %>%
  ggplot() + 
  aes(x = abbr,
      fill = tumor_descriptor) + 
  geom_bar(color = "black", size = 0.25) + 
  scale_fill_manual(values = tumor_descriptor_palette) + 
  facet_wrap(~bh_strip, 
             scales = "free", 
             nrow = 3) +
  labs(
    x = "Cancer group",
    y = "Number of samples", 
    fill = "Tumor descriptor"
  ) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(size = 12, 
                                   angle = 45, 
                                   hjust = 0.8, 
                                   vjust = 0.9),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12))

#Save!
ggsave(filename = file.path(main_output_dir,
                            "tumor_descriptor_proportion_panel.pdf"),
       plot = descriptor_plot,
       width = 10,
       height = 12)


#### Assay types bar plots -----------------------------------------------------

# "Palettes" for the pattern
pattern_key <- c(
  Both = "none",
  `RNA-Seq` = "circle",
  `DNA-Seq` = "stripe"
)

density_key <- c(
  Both = 0.2,
  `RNA-Seq` = 0.7,
  `DNA-Seq` = 0.2
)

spacing_key <- c(
  Both = 0.05,
  `RNA-Seq` = 0.01,
  `DNA-Seq` = 0.05
)


# Create a bar plot that shows the sample size, where the stacked pattern
# shows the experimental strategy
exp_strat_plot <- experimental_strategy_df %>%
  # Order based on display order
  mutate(broad_histology_display = factor(broad_histology_display,
                                          levels = broad_histologies)) %>%
  ggplot(aes(x = cancer_group_display,
             y = n)) +
  geom_bar_pattern(aes(fill = cancer_group_display,
                       pattern = experimental_strategy),
                   color = "#666666",
                   pattern_fill = "#FFFFFF",
                   pattern_color = "#333333",
                   pattern_alpha = 1,
                   stat = "identity") +
  geom_text(aes(y = y_coord,
                label = cancer_group_n)) +
  scale_fill_manual(values = cancer_group_palette) +
  scale_pattern_manual(values =  pattern_key) +
  scale_pattern_density_manual(values = density_key) +
  scale_pattern_spacing_manual(values = spacing_key) +
  theme_pubr() +
  facet_wrap(~ broad_histology_display,
             nrow = 3,
             scales = "free") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 9),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  labs(x = "Cancer Group",
       y = "Sample Size") +
  guides(fill = FALSE,
         color = FALSE,
         pattern = FALSE)

# Save the assay type plot as a supplemental panel
ggsave(filename = file.path(supp_output_dir, "assay_type_stacked_panel.pdf"),
       plot = exp_strat_plot,
       width = 14,
       height = 21)

# We need a little ore control over this legend to make sure the patterns are
# legible
as_ggplot(get_legend(
  experimental_strategy_df %>%
    filter(broad_histology_display == "Diffuse astrocytic and oligodendroglial tumor") %>%
    ggplot(aes(x = cancer_group_display,
               y = n)) +
    geom_bar_pattern(aes(fill = cancer_group_display,
                         pattern = experimental_strategy),
                     color = "#666666",
                     pattern_fill = "#FFFFFF",
                     pattern_color = "#333333",
                     stat = "identity") +
    scale_pattern_manual(values = pattern_key) +
    scale_pattern_density_manual(values = density_key) +
    scale_pattern_spacing_manual(values = spacing_key) +
    theme_classic() +
  labs(pattern = "Assay") +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = FALSE,
         color = FALSE)
))  %>%
  ggsave(filename = file.path(supp_output_dir, "assay_pattern_legend.pdf"),
         height = 2,
         width = 2)
