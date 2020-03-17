# Sample Distribution Figure
#
# 2020
# Chante Bethell for ALSF - CCDL
#
# This script is intended to run steps needed to create Figure 1.

# Load in libraries
library(dplyr)
library(colorspace)
library(scales)
library(treemapify)

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

# Create a treemap of broad histology, short histology, and integrated diagnosis
# Number of palettes needed
n <- length(unique(plots_df$level3))

# Now calculate the colors for each data point
plots_df2 <- plots_df %>%
  mutate(index = as.numeric(factor(plots_df$level3))- 1) %>%
  group_by(index) %>%
  mutate(
    max_size = max(size),
    color = gradient_n_pal(
      sequential_hcl(
        6,
        h = 360 * index[1]/n,
        c = c(45, 20),
        l = c(30, 80),
        power = .5)
    )(size/max_size)
  )

# Plot the treemap
treemap <- ggplot(plots_df2, aes(area = size, fill = color, label=level3, subgroup=level2, subgroup2=level1)) +
    geom_treemap() +
    geom_treemap_subgroup_border(colour="white") +
    geom_treemap_text(fontface = "italic",
                      colour = "white",
                      place = "centre",
                      grow = F,
                      reflow=T) +
    geom_treemap_subgroup_text(place = "top",
                               grow = T,
                               reflow = T,
                               alpha = 0.6,
                               colour = "#FAFAFA",
                               min.size = 0) +
    geom_treemap_subgroup2_text(place = "centre",
                                grow = T,
                                alpha = 0.8,
                                colour = "#FAFAFA",
                                min.size = 0) +
    scale_fill_identity()

## TODO: Re-run Github Contributions plot/table here -- for now we will define
## this plot as NULL
github_contributions_plot <- NULL

## TODO: Re-run or load in plots of the project features and assays -- for now
## we will define these plots as NULL
project_assays_plot <- NULL
project_features_plot <- NULL

#### Assemble multipanel plot -------------------------------------------------
ggpubr::ggarrange(ggpubr::ggarrange(treemap,
                                    project_features_plot,
                                    ncol = 2,
                                    labels = c("A", "C"),
                                    widths = c(2, 1.6),
                                    font.label = list(size = 20)
),
ggpubr::ggarrange(project_assays_plot,
                  github_contributions_plot,
                  ncol = 2,
                  labels = c("B", "D"),
                  widths = c(1.6, 2),
                  font.label = list(size = 20)
),
nrow = 2,
heights = c(2, 2)
)

# Save to PNG
ggplot2::ggsave(file.path(output_dir, "fig1-openpbta-distribution.png"),
                width = 16.5, height = 15,
                units = "in"
)
