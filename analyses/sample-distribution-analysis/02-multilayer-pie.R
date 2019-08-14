# This script creates a treemap and a multilayer pie chart to represent the
# broad histologies, short histologies, and molecular subtypes within the dataset.
#
# This script uses the packages sunburstR, d3r, and treemap to produce the
# visualizations.
#
# Chante Bethell for CCDL 2019

# Load/install packages. This will be removed once a Dockerfile is established.
source("00-install-packages.R", print.eval = FALSE)

# Load source script
source("01-filter-across-types.R")

# Set path to plots directory
outputDir <- file.path("/plots")

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Create a colorblind-friendly color vector
color <- colorblindr::palette_OkabeIto

sun_df <- df2 %>%
  dplyr::select(broad_histology, short_histology, disease_type_new) %>%
  dplyr::filter(!is.na(disease_type_new)) %>%
  dplyr::filter(!is.na(broad_histology)) %>%
  dplyr::filter(!is.na(short_histology)) %>%
  dplyr::group_by(broad_histology)

sun_df_sub <-
  sun_df[, c("broad_histology", "short_histology", "disease_type_new")]
sun_df_sub$nodes <-
  paste0(sun_df_sub[[1]], ",", sun_df_sub[[2]], ",", sun_df_sub[[3]])

# Make node patterns
names(sun_df_sub) <- c("level1", "level2", "level3", "nodes")

# Sum each unique combination
counts <- sun_df_sub %>%
  dplyr::group_by(nodes) %>%
  dplyr::summarise(size = n())

counts <- counts[!grepl("^NA-", counts$nodes),]

# Create final data.frame prepped for treemap and sunburst functions?
final <- merge(sun_df_sub, counts, all.x = T)
col_order <- c("level1", "level2", "level3", "size", "nodes")
final <- final[, col_order]
final$nodes <- NULL

# Save to tsv file
readr::write_tsv(final, file.path("results", "sunburst_plot_df.tsv"))

tm <-
  treemap::treemap(
    final,
    index = c("level1", "level2", "level3"),
    vSize = "size",
    vColor = color,
    draw = TRUE
  )

tmnest <-
  d3r::d3_nest(tm$tm[, c("level1", "level2", "level3", "vSize")],
               value_cols = c("vSize"))

interactive_tm <-
  d3treeR::d3tree(tm,
                  rootname = "Cancer Histologies Treemap",
                  width = 1200,
                  height = 700)

sun_plot <-
  sunburstR::sunburst(
    data = tmnest,
    valueField = "vSize",
    count = TRUE,
    sumNodes = FALSE,
    colors = color
  )

p <- sunburstR::sund2b(tmnest, colors = color, valueField = "vSize")
print(p)

mapview::mapshot(interactive_tm, url = paste0(outputDir, "/histology-treemap.html"))
mapview::mapshot(p, url = paste0(outputDir, "/histology-pie.html"))

sessionInfo()


