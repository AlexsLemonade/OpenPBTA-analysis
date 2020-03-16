# Sample Distribution Figure
#
# 2020
# Chante Bethell for ALSF - CCDL
#
# This script is intended to run steps needed to create Figure 1.

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pngs")

#### Declare relative file paths of modules used for this figure --------------
sample_distribution_dir <- file.path(
  root_dir,
  "analyses",
  "sample-distribution-analysis"
)

## TODO: Define the file path to the directory containing the data for the
##        contributions plot

#### Read in the associated results -------------------------------------------

# Read in disease expression file from `01-filter-across-types.R`
disease_expression <-
  readr::read_tsv(file.path(sample_distribution_dir, "results", "disease_expression.tsv"))

# Read in plots data.frame file from `02-multilaye-plots.R`
plots_df <-
  readr::read_tsv(file.path(sample_distribution_dir, "results", "plots_df.tsv"))

# Reorder the columns to be displayed in descending order by count on the plot
disease_expression$integrated_diagnosis <- with(disease_expression,
                                                reorder(integrated_diagnosis, -count))

#### Re-run the individual plots ----------------------------------------------

# Create a bar plot of sample distribution across cancer types
gg_types_bar <- disease_expression %>%
  ggplot2::ggplot(ggplot2::aes(x = integrated_diagnosis, y = count, fill = count)) +
  ggplot2::geom_col() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Cancer Types", y = "Count",
                title = "Sample Distribution Across Cancer Types") +
  ggplot2::scale_y_continuous(breaks = seq(0, 500, by = 100)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 75,
    hjust = 1,
    size = 8
  ),
  panel.grid = ggplot2::element_blank()) +
  ggplot2::geom_text(nudge_y = 5.5, size = 1,
                     ggplot2::aes(label = paste0(disease_expression$percent)))

# Create a treemap -- if we decide to use the treemap this would need to be
# transformed into a ggplot object for `ggarrange` function in the section below
# tm <-
#   treemap::treemap(
#     plots_df,
#     index = c("level1", "level2", "level3"),
#     vSize = "counter",
#     vColor = color,
#     draw = TRUE
#   )$draw

## TODO: Re-run Github Contributions plot/table here

#### Assemble multipanel plot -------------------------------------------------
ggpubr::ggarrange(ggpubr::ggarrange(gg_types_bar,
                                    ncol = 2,
                                    labels = c("A", ""),
                                    widths = c(2, 1),
                                    font.label = list(size = 22)
),
nrow = 2,
heights = c(2, 2)
)

# Save to PNG
ggplot2::ggsave(file.path(output_dir, "fig1-openpbta-distribution.png"),
                width = 16.5, height = 15,
                units = "in"
)
