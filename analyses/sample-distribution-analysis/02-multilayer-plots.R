# This script creates a treemap and a multilayer pie chart to represent the
# broad histologies, display_group, and molecular subtypes within the dataset.
#
# This script uses the packages sunburstR, d3r, and treemap to produce the
# visualizations.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript analyses/sample-distribution-analysis/02-multilayer-plots.R

# Load in libraries
library(ggplot2)
library(colorspace)
library(scales)
library(treemapify)

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define the file.path to output directories
output_dir <- file.path(root_dir, "analyses", "sample-distribution-analysis")
results_dir <- file.path(output_dir, "results")
plots_dir <- file.path(output_dir, "plots")

# Read in dataset
histologies_df <- readr::read_tsv(file.path(root_dir, "data",
                                            "pbta-histologies.tsv"), guess_max = 10000)

# Read in histology standard color palette for project
histology_label_mapping <- readr::read_tsv(
  file.path(root_dir,
            "figures",
            "palettes", 
            "histology_label_color_table.tsv")) %>% 
  # Select just the columns we will need for plotting
  dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes)

# Create a colorblind-friendly color vector
color <- colorblindr::palette_OkabeIto

# Create final data.frame prepped for treemap and sunburst functions
final_df <- histologies_df %>%
  dplyr::filter(sample_type == "Tumor",
                composition == "Solid Tissue") %>%
  # Join on the color codes
  dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID") %>% 
  # Reorder display_group based on display_order
  dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order)) %>%
  # Get distinct based on participant IDs
  dplyr::distinct(Kids_First_Participant_ID, 
                  broad_histology, 
                  display_group, 
                  harmonized_diagnosis, 
                  hex_codes) %>% 
  # Select our 3 columns of interest
  dplyr::select(broad_histology, display_group, harmonized_diagnosis, hex_codes) %>%
  # Remove any row that has an NA
  dplyr::filter(complete.cases(.)) %>%
  # Group by all 3 columns in order to count
  dplyr::group_by(broad_histology, display_group, harmonized_diagnosis, hex_codes) %>%
  # Add the count to a column named size
  dplyr::add_count(name = "size") %>%
  # Place the value 1 in a column named counter for treemap and sunburt plots
  dplyr::mutate(counter= c(1)) %>%
  # Change the column names
  dplyr::rename(level1 = broad_histology,
                level2 = display_group,
                level3 = harmonized_diagnosis) %>%
  # Reorder the rows according to the 3 levels
  dplyr::arrange(level1, level2, level3) %>%
  # tbl_df -> data.frame
  as.data.frame() 

# Save to tsv file
readr::write_tsv(final_df, file.path(results_dir, "plots_df.tsv"))

# Plot the treemap
treemap <-
  ggplot(
    final_df,
    aes(
      area = size,
      fill = hex_codes,
      label = level2,
      subgroup = level3
    )
  ) +
  geom_treemap() +
  geom_treemap_subgroup_border(colour = "white") +
  geom_treemap_text(
    fontface = "italic",
    colour = "white",
    place = "topleft",
    alpha = 0.4,
    grow = F,
    reflow = T,
    size = 16
  ) +
  geom_treemap_subgroup_text(
    place = "bottomright",
    grow = T,
    reflow = T,
    alpha = 0.6,
    colour = "#FAFAFA",
    min.size = 0
  ) +
  theme(legend.position = "none") +
  scale_fill_identity()

# Save treemap
ggsave(
  treemap,
  file = file.path(plots_dir, "distribution_across_cancer_types_treemap.pdf"),
  width = 22,
  height = 10
)

# Create a treemap (for interactive treemap)
tm <-
  treemap::treemap(
    final_df,
    index = c("level1", "level2", "level3"),
    vSize = "counter",
    vColor = color,
    draw = TRUE
  )$tm

# Convert the tm data.frame into a d3.js hierarchy object which is needed
# for the sund2b plot
tmnest <-
  d3r::d3_nest(tm[, c("level1", "level2", "level3", "vSize")],
               value_cols = c("vSize"))

# Create an interactive treemap
interactive_tm <-
  d3treeR::d3tree(tm,
                  rootname = "Cancer Histologies Treemap",
                  width = 1200,
                  height = 700)

# Create a sunburst plot
sun_plot <-
  sunburstR::sunburst(
    data = tmnest,
    valueField = "vSize",
    count = TRUE,
    sumNodes = FALSE,
    colors = color
  )

# Create an interactive sund2b plot
p <- sunburstR::sund2b(tmnest, colors = color, valueField = "vSize")

# Create HTML outputs for the interactive plots
mapview::mapshot(interactive_tm, url = file.path(plots_dir,
                                                 "histology-treemap.html"))
mapview::mapshot(p, url = file.path(plots_dir, "histology-pie.html"))
