# S. Spielman for ALSF CCDL 2022
#
# Makes a several panels from the `mutational-signatures` module, including:
## Panel 3E, adapted from `analyses/mutational-signatures/07-plot_cns_fit.Rmd`
## Panel S4A, adapted from `analyses/mutational-signatures/07-plot_cns_fit.Rmd`
## Panel S4B, adapted from `analyses/mutational-signatures/07-plot_cns_fit.Rmd`


library(tidyverse)

# Set seed for sina and jitter plot reproducibility
set.seed(2022)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directories
output_dir_fig3 <- file.path(root_dir, "figures", "pdfs", "fig3", "panels")
if (!dir.exists(output_dir_fig3)) {
  dir.create(output_dir_fig3, recursive = TRUE)
}
output_dir_figS4 <- file.path(root_dir, "figures", "pdfs", "supp", "figs4", "panels")
if (!dir.exists(output_dir_figS4)) {
  dir.create(output_dir_figS4, recursive = TRUE)
}

# Zenodo table directory and output file names
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
fig3e_csv <- file.path(zenodo_tables_dir, "figure-3e-data.csv")
figS4a_csv <- file.path(zenodo_tables_dir, "figure-S4a-data.csv")
figS4b_csv <- file.path(zenodo_tables_dir, "figure-S4b-data.csv")


# Directory with result data to input for plot
input_dir <- file.path(root_dir,
                       "analyses",
                       "mutational-signatures",
                       "results")


# Cancer groups we are showing in panel 3E,
#  in same order as 3A but _without_ "Other"
cancer_group_display_order <- c(
  "Diffuse midline glioma",
  "Other high-grade glioma",
  "Pilocytic astrocytoma",
  "Ganglioglioma",
  "Pleomorphic xanthoastrocytoma",
  "Other low-grade glioma",
  "Medulloblastoma",
  "Atypical Teratoid Rhabdoid Tumor",
  "Other embryonal tumor",
  "Ependymoma",
  "Craniopharyngioma",
  "Meningioma"
)


## Input and output files -------------------------------------
meta_file    <- file.path(root_dir, "data", "pbta-histologies.tsv")
palette_file <- file.path(root_dir, "figures", "palettes", "broad_histology_cancer_group_palette.tsv")
tumor_palette_file <- file.path(root_dir, "figures", "palettes", "tumor_descriptor_palette.tsv")
signatures_file  <- file.path(input_dir, "deconstructsigs_exposures_merged.tsv")

sina_iqr_pdf <- file.path(output_dir_fig3, "exposures_sina_IQR.pdf")
barplot_pdf <- file.path(output_dir_figS4, "exposures_per_sample_barplot.pdf")
sig1_descriptor_pdf <- file.path(output_dir_figS4, "signature1_tumor-descriptor_cancer-groups.pdf")



## Prepare data for plots -----------------------
meta       <- read_tsv(meta_file, guess_max = 10000)
palette_df <- read_tsv(palette_file)
tumor_palette_df   <- read_tsv(tumor_palette_file)
signatures_results <- read_tsv(signatures_file) %>%
  # Keep only the groups of interest
  filter(cancer_group_display %in% cancer_group_display_order)


# Arrange cancer_group_display_n in the `cancer_group_display_order` order
# First, find the right order
cancer_group_display_n_order <- signatures_results %>%
  mutate(cancer_group_display = fct_relevel(cancer_group_display,
                                            cancer_group_display_order)) %>%
  arrange(cancer_group_display) %>%
  distinct(cancer_group_display_n) %>%
  pull(cancer_group_display_n)

# Apply the right order
signatures_results <- signatures_results %>%
  mutate(cancer_group_display_n = fct_relevel(cancer_group_display_n,
                                              cancer_group_display_n_order)
  )


# Make panel 3E --------------------------------

exposures_sina_IQR <- ggplot(signatures_results) +
  aes(x = signature,
      y = exposure,
      color = cancer_group_hex) +
  ggforce::geom_sina(size = 0.35) +
  geom_boxplot(outlier.size = 0,
               size = 0.2,
               color = "black",
               alpha = 0,
               # remove whiskers
               coef = 0) +
  facet_wrap(~cancer_group_display_n, nrow = 2) +
  scale_color_identity() +
  labs(
    x = "RefSig signature",
    y = "Signature weights across tumors"
  ) +
  ggpubr::theme_pubr() +
  cowplot::panel_border() +
  theme(
    # angle needed to fit "Other" category
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 5.85),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none"
  )

# Sizing is based on full figure compilation:
ggsave(sina_iqr_pdf,
       exposures_sina_IQR, width = 7.2, height = 3,
       useDingbats = FALSE)

# Export zenodo CSV
signatures_results %>%
  # reorder columns so sample id first
  dplyr::select(Kids_First_Biospecimen_ID, everything()) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # remove \n from cancer_group_display_n so that CSV is properly formatted
  dplyr::mutate(cancer_group_display_n = stringr::str_replace(cancer_group_display_n, "\n", " ")) %>%
  # export
  readr::write_csv(fig3e_csv)


# Make panel S4A --------------------------------

# Make a column for exposure of interest to order samples by
samples_in_order <- signatures_results %>%
  filter(signature == 1) %>%
  select(exposure, Kids_First_Biospecimen_ID) %>%
  arrange(-exposure) %>%
  pull(Kids_First_Biospecimen_ID)

signatures_results_fct <- signatures_results %>%
  mutate(Kids_First_Biospecimen_ID = factor(Kids_First_Biospecimen_ID,
                                            levels = samples_in_order))
sample_exposures_barplot <- ggplot(signatures_results_fct) +
  aes(x = Kids_First_Biospecimen_ID,
      y = exposure,
      fill = signature) +
  geom_col(color = "black", size = 0.05) +
  labs(x = "Tumor",
       y = "Signature weight",
       fill = "RefSig Signature") +
  # Only 1:8 palette. The "Other" gets filled as white by including `color="black"` in the geom
  # There will be a warning about palette size and the warning is OK because ^^
  colorblindr::scale_fill_OkabeIto() +
  facet_wrap(~cancer_group_display_n, scales = "free_x", nrow = 3) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title   = element_text(size = 10),
    axis.text.y   = element_text(size = 8),
    strip.text   = element_text(size = 8.5),
    axis.line = element_line(size = rel(0.5)),
    axis.ticks.y = element_line(size = rel(0.5)),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

ggsave(barplot_pdf,
       sample_exposures_barplot,
       width = 8, height = 6.5)

# Export zenodo CSV
signatures_results_fct %>%
  # reorder columns so sample id first
  dplyr::select(Kids_First_Biospecimen_ID, everything()) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # remove \n from cancer_group_display_n so that CSV is properly formatted
  dplyr::mutate(cancer_group_display_n = stringr::str_replace(cancer_group_display_n, "\n", " ")) %>%
  # export
  readr::write_csv(figS4a_csv)



## Make panel S4B ---------------------------


# Grab the tumor descriptor colors
tumor_colors <- tumor_palette_df$hex_codes
names(tumor_colors) <- tumor_palette_df$color_names


sig1_descriptor_plot_df <- signatures_results %>%
  filter(signature == 1) %>%
  rename(sig1_proportion = exposure)

sig1_descriptor_plot <- ggplot(sig1_descriptor_plot_df) +
  aes(x = tumor_descriptor,
      y = sig1_proportion,
      color = tumor_descriptor) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  # light guiding line representing 0 exposure
  geom_hline(yintercept = 0, size = 0.15) +
  scale_color_manual(values = tumor_colors) +
  # add in mean +/- SE pointrange
  stat_summary(color = "black", size = 0.3) +
  facet_wrap(~cancer_group_display_n,
             nrow = 2) +
  labs(
    x = "Tumor descriptor",
    y = "Signature 1 Weight"
  ) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.25),
    axis.text.y = element_text(size = 7.5),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 6),
    axis.line = element_line(size = rel(0.5)),
    axis.ticks = element_line(size = rel(0.5)),
    legend.position = "none"
  )


ggsave(sig1_descriptor_pdf,
       sig1_descriptor_plot,
       width = 8, height = 4,
       useDingbats = FALSE)


# Export zenodo CSV
sig1_descriptor_plot_df %>%
  # reorder columns so sample id first
  dplyr::select(Kids_First_Biospecimen_ID, everything()) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # remove \n from cancer_group_display_n so that CSV is properly formatted
  dplyr::mutate(cancer_group_display_n = stringr::str_replace(cancer_group_display_n, "\n", " ")) %>%
  # export
  readr::write_csv(figS4b_csv)

