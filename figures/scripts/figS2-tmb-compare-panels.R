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

# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig2h_csv <- file.path(zenodo_tables_dir, "figure-S2h-data.csv") # PBTA
fig2i_csv <- file.path(zenodo_tables_dir, "figure-S2i-data.csv") # TCGA

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

# Read in manifest to get Primary_diagnosis
tmb_manifest <- read_tsv(
  file.path(data_dir, "pbta-tcga-manifest.tsv")
)

# Join diagnosis information with tmb_tcga
tmb_tcga <- inner_join(tmb_tcga, 
                       select(tmb_manifest, 
                              Tumor_Sample_Barcode = tumorID,
                              Primary_diagnosis))


# Define functions ----------------------------------------------------

# Function to calculate medians and ranks
prepare_data_for_plot <- function(df, grouping_variable = NULL, min_samples = 5) {
  df %>%
    # Group by specified column
    group_by({{grouping_variable}}) %>%
    # Only keep groups with the specified minimum number of samples
    filter(n() > min_samples) %>%
    # Calculate group median
    mutate(
      group_median = median(tmb, na.rm = TRUE),
      group_rank = rank(tmb, ties.method = "first") / n(),
      sample_size = paste0("n = ", n())
    ) %>%
    ungroup() 
}



# Function to plot shared layers of PBTA and TCGA plots
# Note this isn't the entire plot due to tidyeval challenges with facet variables
baseline_plot <- function(df, tcga_color = NULL) {
  
  # Start plot
  p <-  ggplot(df) + 
    aes(
      x = group_rank,
      y = tmb
    ) 
  
  # Add point with or without color 
  if (!(is.null(tcga_color))) {
    p <- p + geom_point(color = tcga_color, alpha = 0.7)
  } else {
    p <- p + geom_point(aes(color = cancer_group_hex), alpha = 0.7)
  }
  
  # Rest of plot
  p + 
    # Add summary line for median
    geom_segment(
      x = 0, xend = 1, color = "black",
      aes(y = group_median, yend = group_median)
    ) +
    xlim(-1.2, 1.2) +
    scale_y_continuous(
      trans = "log1p",
      limits = c(0, 400),
      breaks = c(0, 3, 10, 30, 100, 300)
    ) +
    labs(
      x = "", # "Cancer group",
      y = "Coding mutations per Mb"
    ) +
    ggpubr::theme_pubr() +
    theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 9, angle = 90, hjust = 1),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      legend.position = "none"
    ) 
}



# Prepare data for plotting --------------------------------------------------

# Combine tmb_pbta data with metadata and cancer group information
tmb_pbta_plot_df <- tmb_pbta %>%
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
  # Perform calculations needed for plot
  prepare_data_for_plot(grouping_variable = cancer_group_display) %>%
  # remove "Other" cancer group
  filter(cancer_group_display != "Other") %>%
  # Order cancer groups by median TMB
  mutate(cancer_group_display = str_wrap(cancer_group_display, 18),
         cancer_group_display = fct_reorder(cancer_group_display, tmb, .fun = median)
  ) 

# Prepare tcga data
tmb_tcga_plot_df <- tmb_tcga %>%
  prepare_data_for_plot(grouping_variable = Primary_diagnosis) %>%
  # Order Primary_diagnosis by median TMB
  mutate(
    Primary_diagnosis = str_wrap(Primary_diagnosis, 10),
    Primary_diagnosis = fct_reorder(Primary_diagnosis, tmb, .fun = median),
  )

# Plot the data ----------------------------------------------------------
tmb_pbta_plot <- baseline_plot(tmb_pbta_plot_df) +
  facet_wrap(~cancer_group_display + sample_size, nrow = 1, strip.position = "bottom")  +
  scale_color_identity()
  

tmb_tcga_plot <- baseline_plot(tmb_tcga_plot_df, tcga_color = "gray60") +
  facet_wrap(~Primary_diagnosis + sample_size, nrow = 1, strip.position = "bottom") 

# Export both plots ---------------------------------
ggsave(pbta_tmb_cdf_pdf, tmb_pbta_plot, width = 10, height = 6, 
useDingbats = FALSE)
ggsave(tcga_tmb_cdf_pdf, tmb_tcga_plot, width = 5, height = 6.15, # because labels are diff sizes, making this 6.25 matches the other plot at 6
useDingbats = FALSE)  


# Export Zenodo CSV files --------------------

# S2H
tmb_pbta_plot_df %>%
  # reorder columns and remove unnecessary columns, but keep hex since identity was plotted
  dplyr::select(Kids_First_Biospecimen_ID, everything(),
                -cancer_group, -cancer_group_abbreviation) %>%
  # remove \n from wrapped display column
  dplyr::mutate(cancer_group_display = stringr::str_replace_all(cancer_group_display, "\n", " ")) %>%
  # remove "n = " from sample size column
  dplyr::mutate(sample_size = stringr::str_replace(sample_size, "n = ", "")) %>%
  # arrange
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # export
  readr::write_csv(fig2h_csv)


# S2I
tmb_tcga_plot_df %>%
  # reorder columns
  dplyr::select(Tumor_Sample_Barcode, Primary_diagnosis, everything()) %>%
  # remove \n from wrapped display column
  dplyr::mutate(Primary_diagnosis = stringr::str_replace(Primary_diagnosis, "\n", " ")) %>%
  # remove "n = " from sample size column
  dplyr::mutate(sample_size = stringr::str_replace(sample_size, "n = ", "")) %>%
  # arrange on Primary_diagnosis in this case since it is not PBTA data
  dplyr::arrange(Primary_diagnosis) %>%
  # export
  readr::write_csv(fig2i_csv)
    



