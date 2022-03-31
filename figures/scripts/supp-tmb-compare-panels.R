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

# Read in clinical data
metadata <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                     guess_max = 10000)

# Read in histologies palette file
histologies_palette_df <- read_tsv(file.path(root_dir, "figures", 
                                             "palettes", "broad_histology_cancer_group_palette.tsv"))


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
    as_tibble() %>%
    select(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode,
           tmb) %>%
    inner_join(
      select(metadata,
             Kids_First_Biospecimen_ID,
             cancer_group)
      ) %>%
    inner_join(
      select(histologies_palette_df, 
             contains("cancer_group")
      )
    ) %>%
    drop_na(cancer_group_display) %>%
    # Group by specified column
    group_by(cancer_group_display) %>%
    # Only keep groups with the specified minimum number of samples
    filter(n() > min_samples) %>%
    # Calculate group median
    mutate(
      cancer_group_median = median(tmb, na.rm = TRUE),
      cancer_group_rank = rank(tmb, ties.method = "first") / n(),
      sample_size = paste0("n = ", n())
    ) %>%
    ungroup() %>%
    # Order cancer groups in frequency, but ensure Other is still last
    mutate(cancer_group_display = str_wrap(cancer_group_display, 18),
           cancer_group_display = fct_infreq(cancer_group_display),
           # reverse since infreq does most --> least
           cancer_group_display = fct_rev(cancer_group_display),
           # move "Other" to the end
           cancer_group_display = fct_relevel(cancer_group_display, "Other", after = Inf)
    ) 
}




plot_tmb <- function(df, ylim, ybreaks) {
  ggplot(df) +
    aes(
      x = cancer_group_rank,
      y = tmb,
      color = cancer_group_hex
    ) +
    geom_point() +
    # Add summary line for median
    geom_segment(
      x = 0, xend = 1, color = "black",
      aes(y = cancer_group_median, yend = cancer_group_median)
    ) +
    scale_color_identity() +
    facet_wrap(~ cancer_group_display + sample_size, nrow = 1, strip.position = "bottom") +
    labs(
      x = "Cancer group",
      y = "Coding mutations per Mb"
    ) +
    # Transform to log10 make non-log y-axis labels
    scale_y_continuous(
      trans = "log1p",
      limits = ylim,
      breaks = ybreaks
    ) +
    xlim(-1.2, 1.2) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 11, angle = 90, hjust = 1),
      strip.background = ggplot2::element_rect(fill = NA, color = NA)
    ) 
}

# Run preparation function
tmb_pbta_plot_df <- prepare_data_for_plot(tmb_pbta)
tmb_tcga_plot_df <- prepare_data_for_plot(tmb_tcga) # BUG HERE


plot_tmb(tmb_pbta_plot_df,
         c(0, 400),
         c(0, 3, 10, 30, 100, 300))


plot_tmb(tmb_tcga_plot_df,
         c(0, 400),
         c(0, 3, 10, 30, 100, 300))



# Plot ------------------
tmb_pbta_plot <- ggplot(tmb_pbta_plot_df) +
  aes(
    x = cancer_group_rank,
    y = tmb,
    color = cancer_group_hex
  ) +
  geom_point() +
  # Add summary line for median
  geom_segment(
    x = 0, xend = 1, color = "black",
    aes(y = cancer_group_median, yend = cancer_group_median)
  ) +
  scale_color_identity() +
  facet_wrap(~ cancer_group_display + sample_size, nrow = 1, strip.position = "bottom") +
  labs(
    x = "Cancer group",
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
    strip.text = ggplot2::element_text(size = 11, angle = 90, hjust = 1),
    strip.background = ggplot2::element_rect(fill = NA, color = NA)
  ) 


# Export both plots
ggsave(pbta_tmb_cdf_pdf, pbta_plot, width = 22, height = 4)
ggsave(tcga_tmb_cdf_pdf, tcga_plot, width = 16, height = 4)







