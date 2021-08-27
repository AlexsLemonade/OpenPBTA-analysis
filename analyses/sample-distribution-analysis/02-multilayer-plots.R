# This script creates a treemap and a multilayer pie chart to represent the
# broad histologies, cancer_group, and molecular subtypes within the dataset.
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
palette_dir <- file.path(root_dir, "figures", "palettes")

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
  dplyr::select(Kids_First_Biospecimen_ID, cancer_group_order, cancer_group_hex_codes)

# Use histology_label_color_table.tsv for broad_histology and cancer_group hex codes
histologies_color_key_df <- readr::read_tsv(file.path(palette_dir,"histology_label_color_table.tsv"),
                                            col_types = readr::cols())

# Create final data.frame prepped for treemap and sunburst functions
final_df <- histologies_df %>%
  dplyr::filter(sample_type == "Tumor",
                composition == "Solid Tissue") %>%
  # Join on the color codes
  dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID") %>% 
  # Reorder cancer_group based on cancer_group_order
  dplyr::mutate(cancer_group = forcats::fct_reorder(cancer_group, cancer_group_order)) %>%
  # Extract WGS, WXS, RNA-Seq
  # To-do targeted sequencing?
  reshape2::dcast(Kids_First_Participant_ID + 
                    broad_histology +
                    cancer_group + 
                    tumor_descriptor + 
                    germline_sex_estimate +
                    cancer_group_hex_codes ~ experimental_strategy,
                  fun.aggregate = function(x){as.integer(length(x)>0)}) %>%
  dplyr::mutate(WGS = dplyr::if_else(WGS==1,"WGS","Not Available"),
                WXS = dplyr::if_else(WXS==1,"WXS","Not Available"),
                `RNA-Seq`= dplyr::if_else(`RNA-Seq`==1,"RNA-Seq","Not Available")) %>%
  # Get distinct based on participant IDs
  dplyr::distinct(Kids_First_Participant_ID, 
                  WGS,
                  WXS,
                  `RNA-Seq`,
                  broad_histology,
                  cancer_group, 
                  tumor_descriptor,
                  germline_sex_estimate, 
                  cancer_group_hex_codes) %>% 
  # Select our columns of interest
  dplyr::select(broad_histology,
                cancer_group,
                WGS,
                WXS, 
                `RNA-Seq`, 
                tumor_descriptor,
                germline_sex_estimate,
                cancer_group_hex_codes) %>%
  # Remove any row that has an NA
  dplyr::filter(complete.cases(.)) %>%
  # Group by all columns in order to count
  dplyr::group_by( broad_histology,
                   cancer_group,
                   WGS,
                   WXS,
                   `RNA-Seq`,
                   tumor_descriptor,
                   germline_sex_estimate,
                   cancer_group_hex_codes) %>%
  # Add the count to a column named size
  dplyr::add_count(name = "size") %>%
  # Place the value 1 in a column named counter for treemap and sunburt plots
  dplyr::mutate(counter= c(1)) %>%
  # Change the column names
  dplyr::rename(level1 = broad_histology,
                level2 = cancer_group,
                level3 = WGS,
                level4 = WXS,
                level5 = `RNA-Seq`,
                level6 = tumor_descriptor,
                level7 = germline_sex_estimate) %>%
  # Reorder the rows according to the 3 levels
  dplyr::arrange(level1, level2, level3, level4, level5, level6, level7) %>%
  # tbl_df -> data.frame
  as.data.frame() 

# Save to tsv file
readr::write_tsv(final_df, file.path(results_dir, "plot_df.tsv"))

# Create a treemap (for interactive treemap)
tm <-
  treemap::treemap(
    final_df,
    index = c("level1", "level2", "level3",
              "level4", "level5", "level6"),
    vSize = "counter",
    draw = TRUE
  )$tm

# Update colors
level1 <- subset(tm, level == 1)
# merge to get hex_codes 
level1 <- level1 %>%
  dplyr::select(-color)%>%
  dplyr::left_join(histologies_color_key_df
                          , by= c( "level1"="broad_histology")) %>%
  dplyr::rename(color=hex_codes)%>%
  dplyr::select(level1, color)
level2 <- subset(tm, level == 2)
# merge to get hex_codes
level2 <- level2 %>%
  dplyr::select(-color)%>%
  dplyr::left_join(histologies_color_key_df
                   , by= c( "level2"="cancer_group")) %>%
  dplyr::rename(color=cancer_group_hex_codes)%>%
  dplyr::select(level2, color)
level3 <- subset(tm, level == 3)
level3$color <- ifelse(level3$level3 == "WGS", "#E5E5E5", "#FFFFFF")
level4 <- subset(tm, level == 4)
level4$color <- ifelse(level4$level4 == "WXS", "#B3B3B3", "#FFFFFF")
level5 <- subset(tm, level == 5)
level5$color <- ifelse(level5$level5 == "RNA-Seq", "#CCCCCC", "#FFFFFF")
level6 <- subset(tm, level == 6)
level7 <- subset(tm, level == 7)

new.tm <- dplyr::bind_rows(list(level1, level2, level3, level4, level5, level6, level7))

# Convert the tm data.frame into a d3.js hierarchy object which is needed
# for the sund2b plot
tmnest <-
  d3r::d3_nest(tm[, c("level1", "level2", "level3",
                      "level4", "level5", "level6",
                      "vSize")],
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
    sumNodes = FALSE
  )

# Create an interactive sund2b plot
p <- sunburstR::sund2b(tmnest, valueField = "vSize")

# Create HTML outputs for the interactive plots
mapview::mapshot(p, url = file.path(plots_dir, "histology-pie.html"))
