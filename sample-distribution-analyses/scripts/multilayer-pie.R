# This script creates a treemap and a multilayer pie chart to represent the broad histologies, short histologies, and molecular subtypes within the dataset. 
# This script uses the packages sunburstR, d3r, and treemap to produce the visualizations. 
#
# CCDL - 2019 
# Chante Bethell 

# Load/install packages. This will be removed once a Dockerfile is established. 
source(file.path("..", "scripts", "install-packages.R"), print.eval = FALSE)

#Load source script
source(file.path("pbta-analysis-5.R"))

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Create a colorblind-friendly color vector 
color <- colorblindr::palette_OkabeIto

sun.df <- df2 %>% 
  select(broad_histology, short_histology, disease_type_new) %>%
  filter(!is.na(disease_type_new)) %>%
  filter(!is.na(broad_histology)) %>%
  filter(!is.na(short_histology)) %>%
  group_by(broad_histology)

sun.df.sub <- sun.df[,c("broad_histology", "short_histology", "disease_type_new")]
sun.df.sub$nodes <- paste0(sun.df.sub[[1]], ",", sun.df.sub[[2]], ",", sun.df.sub[[3]])

# Make node patterns
names(sun.df.sub) <- c("level1", "level2", "level3", "nodes")

# Sum each unique parent-node-leaf combination
counts <- sun.df.sub %>% 
  group_by(nodes) %>%
  dplyr::summarise(size = n())
# Remove with no assay
counts <- counts[!grepl("^NA-", counts$nodes), ]

final <- merge(sun.df.sub, counts, all.x = T)
col.order <- c("level1", "level2", "level3",
               "size", "nodes")
final <- final[, col.order]
final$nodes <- NULL 

# Save to tsv file
write_tsv(final, file.path("..", "results", "sunburst_plot_df.tsv"))

tm <- treemap(final, index=c("level1", "level2", "level3"), 
              vSize="size", vColor= color, draw=TRUE)
tmnest <- d3_nest(tm$tm[, c("level1", "level2", "level3", "vSize")],
                  value_cols = c("vSize"))
sun.plot <- sunburst(
  data = tmnest,
  valueField = "vSize",
  count = TRUE,
  sumNodes = FALSE,
  colors = color
)

p <- sund2b(
  tmnest,
  colors = color,
  valueField = "vSize"
)

p


sessionInfo()


