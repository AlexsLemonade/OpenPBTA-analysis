# Sample Distribution Figure
#
# 2020
# Chante Bethell for ALSF - CCDL
#
# This script is intended to run steps needed to create Figure 1.

# Load in libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(scales)
library(treemapify)
library(patchwork)

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

# Read in plots data.frame file from `02-multilayer-plots.R`
plots_df <-
  readr::read_tsv(file.path(sample_distribution_dir, "results", "plots_df.tsv")) %>% 
  dplyr::select(-hex_codes)

# Reorder the columns to be displayed in descending order by count on the plot
disease_expression$harmonized_diagnosis <- with(disease_expression,
                                                reorder(harmonized_diagnosis, -count))

# Read in the histology color palette
histology_label_mapping <- readr::read_tsv(
  file.path(root_dir,
            "figures",
            "palettes",
            "histology_label_color_table.tsv")) %>%
  # Select just the columns we will need for making sure the hex_codes are up to date
  dplyr::select(display_group, display_order, hex_codes) %>% 
  dplyr::distinct()

#### Re-run the individual plots ----------------------------------------------

# Create a treemap of broad histology, short histology, and harmonized diagnosis

# Join the color palette for the colors for each short histology value --
# palette is generated in `figures/scripts/color_palettes.R`
plots_df2 <- plots_df %>%
  left_join(histology_label_mapping, by = c("level2" = "display_group")) %>%
  distinct() # Remove the redundant rows from prep for the `treemap` function

# Plot the treemap where level1 is `broad_histology`,
# level2 is `display_group`, and level3 is `harmonized_diagnosis`
treemap <-
  ggplot(
    plots_df2,
    aes(
      area = size,
      fill = hex_codes,
      label = level3,
      subgroup = level1
    )
  ) +
  geom_treemap() +
  geom_treemap_subgroup_border(colour = "white") +
  geom_treemap_text(
    fontface = "italic",
    colour = "white",
    place = "topright",
    alpha = 0.3,
    grow = F,
    reflow = T,
    min.size = 0,
    size = 6
  ) +
  geom_treemap_subgroup_text(
    place = "bottomleft",
    grow = F,
    reflow = T,
    alpha = 0.6,
    colour = "#FAFAFA",
    size = 10
  ) +
  theme(legend.position = "none") +
  scale_fill_identity()

## TODO: Re-run Github Contributions plot/table here -- for now we will define
## this plot as NULL
github_contributions_plot <- NULL

## TODO: Re-run or load in plots of the project features and assays -- for now
## we will define these plots as NULL
project_assays_plot <- NULL
project_features_plot <- NULL

#### Assemble multipanel plot -------------------------------------------------

# Combine plots with patchwork
# Layout of the four plots will be two over the other two
# (2 columns and 2 rows)
combined_plot <- treemap + project_features_plot +
  project_assays_plot + github_contributions_plot +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = 'A') &
  theme(# add uniform labels
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9))

# Save to PNG
ggplot2::ggsave(file.path(output_dir, "fig1-openpbta-distribution.png"),
                width = 12, height = 8,
                units = "in"
)
