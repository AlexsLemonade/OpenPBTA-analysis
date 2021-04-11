# Mutational Landscape Figure
#
# 2020
# C. Savonen for ALSF - CCDL
#
# Purpose  Run steps needed to create mutational-landscape Figure.
#
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

######## Declare relative file paths of modules used for this figure ###########
snv_callers_dir <- file.path(
  root_dir,
  "analyses",
  "snv-callers"
)
mut_sig_dir <- file.path(
  root_dir,
  "analyses",
  "mutational-signatures"
)
tmb_dir <- file.path(
  root_dir,
  "analyses",
  "tmb-compare"
)

# Import specialized functions from mutational-signatures
source(file.path(mut_sig_dir, "util", "mut_sig_functions.R"))
source(file.path(tmb_dir, "util", "cdf-plot-function.R"))

# Read in the color palette
gradient_col_palette <- readr::read_tsv(
  file.path(root_dir, "figures", "palettes", "gradient_color_palette.tsv")
  )

# Won't need NA color this time.
gradient_col_palette <- gradient_col_palette %>%
  dplyr::filter(color_names != "na_color")

###################### Read in associated results ##############################
# Read in PBTA TMB results
tmb_pbta <- data.table::fread(file.path(
  snv_callers_dir,
  "results",
  "consensus",
  "pbta-snv-mutation-tmb-coding.tsv"
)) %>%
  # Remove Panel samples
  dplyr::filter(experimental_strategy != "Panel")

# Read in TCGA TMB results
tmb_tcga <- data.table::fread(file.path(
  snv_callers_dir,
  "results",
  "consensus",
  "tcga-snv-mutation-tmb-coding.tsv"
))

# Read in cosmic signature results
cosmic_sigs_df <- readr::read_tsv(file.path(
  mut_sig_dir,
  "results",
  "cosmic_signatures_results.tsv"
)) %>%
  # Get rid of `.` in signature names and factor in order
  dplyr::mutate(
    signature = gsub("\\.", " ", signature),
    signature = factor(signature, levels = unique(signature))
  )

# Read in nature signature results
nature_sigs_df <- readr::read_tsv(file.path(
  mut_sig_dir,
  "results",
  "nature_signatures_results.tsv"
)) %>%
  # Get rid of `.` in signature names and factor in order
  dplyr::mutate(
    signature = gsub("\\.", " ", signature),
    signature = factor(signature, levels = unique(signature))
  )

###################### Re-run the individual plots #############################
# Make PBTA TMB plot
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
  breaks = c(0, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)
) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 12),
    plot.margin = ggplot2::unit(c(1, 1, -2, 1), "cm")
  )

# Make TCGA plot
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
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_text(size = 11),
    plot.margin = ggplot2::unit(c(1, 1, -4.4, 0), "cm")
  )

# Create cosmic signatures bubble plot
mut_sig_plot_cosmic <- bubble_matrix_plot(cosmic_sigs_df,
  label = "COSMIC",
  color_palette = gradient_col_palette$hex_codes
) +
  ggplot2::theme(
    legend.position = "none",
    axis.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::unit(c(.5, 1, .5, .5), "cm")
  )

# Create nature signatures bubble plot
mut_sig_plot_nature <- bubble_matrix_plot(nature_sigs_df,
  label = "Alexandrov et al, 2013",
  color_palette = gradient_col_palette$hex_codes
) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::unit(c(.5, .5, .5, .5), "cm")
  )

########################### Assemble multipanel plot ###########################
ggpubr::ggarrange(ggpubr::ggarrange(pbta_plot,
  tcga_plot,
  ncol = 2,
  labels = c("A", ""),
  widths = c(2, 1),
  font.label = list(size = 22)
),
ggpubr::ggarrange(mut_sig_plot_cosmic,
  mut_sig_plot_nature,
  ncol = 2,
  labels = c("B", "C"),
  widths = c(1.6, 2),
  font.label = list(size = 22)
),
nrow = 2,
heights = c(1.5, 2)
)

# Save to PNG
ggplot2::ggsave(file.path(output_dir, "fig2-mutational-landscapes.png"),
  width = 16.5, height = 15,
  units = "in"
)
