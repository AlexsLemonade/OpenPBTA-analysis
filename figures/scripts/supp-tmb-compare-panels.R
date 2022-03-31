# S. Spielman for CCDL 2022
#
# Makes pdf panels for supplementary Figure S2, specifically those that are derived from the `tmb-compare` analysis module.

library(tidyverse)


# Directories -------------------------------------------------------------------
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "supp", "figs2", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory
data_dir <- file.path(root_dir, "data")


## Define final PDF files ------------------------------------------------------
pbta_tmb_cdf_pdf <- file.path(output_dir, "pbta_tmb_cdf_plot.pdf")
tcga_tmb_cdf_pdf <- file.path(output_dir, "tcga_tmb_cdf_plot.pdf")

# Read in consensus data from pbta and tcga --------------------------------------------------------------
tmb_pbta <- data.table::fread(file.path(
  data_dir,
  "pbta-snv-consensus-mutation-tmb-coding.tsv"
  )) %>%
  # apparently, this variable is weird when binding but we don't need it for the plot so we'll just remove it.
  select(-region_size) %>%
  filter(experimental_strategy != "Panel")

tmb_tcga <- data.table::fread(file.path(
  data_dir,
  "tcga-snv-mutation-tmb-coding.tsv"
  )) %>%
  select(-region_size)

# Prepare data for plotting ----------------------------------------------------
prepare_data_for_plot <- function(df, min_samples = 5) {
  df %>%
    # We only really need these two variables
    transmute(
      short_histology = str_to_title(short_histology),
      tmb = as.numeric(tmb)
    ) %>%
    # Group by specified column
    group_by(short_histology) %>%
    # Only keep groups with the specified minimum number of samples
    filter(n() > min_samples) %>%
    # Calculate group (short_histology) median
    dplyr::mutate(
      group_median = median(tmb, na.rm = TRUE),
      group_rank = rank(tmb, ties.method = "first") / n(),
      sample_size = paste0("n = ", n())
    ) %>%
    ungroup() %>%
    mutate(short_histology = fct_reorder(short_histology, group_median))
}

# Run preparation function
tmb_pbta_plot_df <- prepare_data_for_plot(tmb_pbta)
tmb_tcga_plot_df <- prepare_data_for_plot(tmb_tcga)

# Plot ------------------
ggplot(tmb_pbta_plot_df) +
  aes(
    x = group_rank,
    y = number
  ) +
  geom_point(color = "blue") +
  # Add summary line for median
  geom_segment(
    x = 0, xend = 1, color = "grey",
    aes(y = group_median, yend = group_median)
  ) +
  # Separate by histology
  facet_wrap(~ group + sample_size, nrow = 1, strip.position = "bottom") +
  labs(
    x = "Short histology",
    y = "Coding mutations per Mb"
  ) +
  # Transform to log10 make non-log y-axis labels
  scale_y_continuous(
    trans = "log1p",
    limits = c(0, 400),
    breaks = c(0, 3, 10, 30, 100, 300)
  ) +
  xlim(-1.2, 1.2) +
  ggpubr::theme_pubr() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.text = ggplot2::element_text(size = 10, angle = 90, hjust = 1),
    strip.background = ggplot2::element_rect(fill = NA, color = NA)
  ) 



# Plot a CDF plot for each of pbta and tcga
pbta_plot <- cdf_plot(
  df = tmb_pbta,
  plot_title = "PBTA",
  num_col = "tmb",
  group_col = "short_histology",
  color = "#3BC8A2",
  n_group = 5,
  x_lim = c(-1.2, 1.2),
  y_lim = c(0, 400),
  x_lab = "",
  y_lab = "Coding Mutations per Mb", 
  breaks = c(0, 3, 10, 30, 100, 300)
) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x =  element_text(size = 7),
    strip.text.x = element_text(size = 8)
  )

tcga_plot <- cdf_plot(
  df = tmb_tcga,
  plot_title = "TCGA (Adult)",
  num_col = "tmb",
  group_col = "short_histology",
  color = "#630882",
  n_group = 5,
  x_lim = c(-1.2, 1.2),
  y_lim = c(0, 400),
  x_lab = "",
  y_lab = "Coding Mutations per Mb",
  breaks = c()
) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x =  element_text(size = 7),
    strip.text.x = element_text(size = 9)
  )

# Export both plots
ggsave(pbta_tmb_cdf_pdf, pbta_plot, width = 22, height = 4)
ggsave(tcga_tmb_cdf_pdf, tcga_plot, width = 16, height = 4)







