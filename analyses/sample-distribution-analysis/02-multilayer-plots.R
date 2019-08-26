# This script creates a treemap and a multilayer pie chart to represent the
# broad histologies, short histologies, and molecular subtypes within the dataset.
#
# This script uses the packages sunburstR, d3r, and treemap to produce the
# visualizations.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory of the repository as follows:
# Rscript analyses/sample-distribution-analysis/02-multilayer-plots.R

# Set paths to output directories
outputDir <- file.path("analyses", "sample-distribution-analysis")
results_dir <- file.path(outputDir, "results")
plots_dir <- file.path(outputDir, "plots")

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Read in dataset 
df2 <- readr::read_tsv(file.path("data",
                                 "pbta-histologies.tsv"))

# Create a colorblind-friendly color vector
color <- colorblindr::palette_OkabeIto

# Create final data.frame prepped for treemap and sunburst functions
final_df <- df2 %>%
  # Select our 3 columns of interest
  dplyr::select(broad_histology, short_histology, disease_type_new) %>%
  # Remove any row that has an NA
  dplyr::filter(complete.cases(.)) %>%
  # Group by all 3 columns in order to count
  dplyr::group_by(broad_histology, short_histology, disease_type_new) %>%
  # Add the count to a column named size
  dplyr::add_count(name = "size") %>%
  # Place the value 1 in a column named counter for treemap and sunburt plots
  dplyr::mutate(counter= c(1)) %>%
  # Change the column names
  dplyr::rename(level1 = broad_histology, 
                level2 = short_histology,
                level3 = disease_type_new) %>%
  # Reorder the rows according to the 3 levels
  dplyr::arrange(level1, level2, level3) %>%
  # tbl_df -> data.frame
  as.data.frame()

# Save to tsv file
readr::write_tsv(final_df, file.path(results_dir, "plots_df.tsv"))

# Create a treemap 
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
mapview::mapshot(interactive_tm, url = file.path(plots_dir, "histology-treemap.html"))
mapview::mapshot(p, url = file.path(plots_dir, "histology-pie.html"))



