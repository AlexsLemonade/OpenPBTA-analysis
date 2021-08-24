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

# Create a colorblind-friendly color vector
color <- colorblindr::palette_OkabeIto

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
                    cancer_group + 
                    tumor_descriptor + 
                    germline_sex_estimate +
                    cancer_group_hex_codes ~ experimental_strategy,
                  fun.aggregate = function(x){as.integer(length(x)>0)}) %>%
  dplyr::mutate(WGS = dplyr::if_else(WGS==1,"Available","Not Available"),
                WXS = dplyr::if_else(WXS==1,"Available","Not Available"),
                `RNA-Seq`= dplyr::if_else(`RNA-Seq`==1,"Available","Not Available")) %>%
  # Get distinct based on participant IDs
  dplyr::distinct(Kids_First_Participant_ID, 
                  WGS,
                  WXS,
                  `RNA-Seq`,
                  cancer_group, 
                  tumor_descriptor,
                  germline_sex_estimate, 
                  cancer_group_hex_codes) %>% 
  # Select our columns of interest
  dplyr::select(cancer_group, WGS, WXS, `RNA-Seq`, tumor_descriptor, germline_sex_estimate, cancer_group_hex_codes) %>%
  # Remove any row that has an NA
  dplyr::filter(complete.cases(.)) %>%
  # Group by all columns in order to count
  dplyr::group_by( cancer_group, WGS, WXS, `RNA-Seq`, tumor_descriptor, germline_sex_estimate, cancer_group_hex_codes) %>%
  # Add the count to a column named size
  dplyr::add_count(name = "size") %>%
  # Place the value 1 in a column named counter for treemap and sunburt plots
  dplyr::mutate(counter= c(1)) %>%
  # Change the column names
  dplyr::rename(level1 = cancer_group,
                level2 = WGS,
                level3 = WXS,
                level4 = `RNA-Seq`,
                level5 = tumor_descriptor,
                level6 = germline_sex_estimate) %>%
  # Reorder the rows according to the 3 levels
  dplyr::arrange(level1, level2, level3, level4, level5, level6) %>%
  # tbl_df -> data.frame
  as.data.frame() 

# Save to tsv file
readr::write_tsv(final_df, file.path(results_dir, "sample_dist_plot_df.tsv"))

# Create a treemap (for interactive treemap)
tm <-
  treemap::treemap(
    final_df,
    index = c("level1", "level2", "level3",
              "level4", "level5", "level6"),
    vSize = "counter",
    vColor = color,
    draw = TRUE
  )$tm

# Convert the tm data.frame into a d3.js hierarchy object which is needed
# for the sund2b plot
tmnest <-
  d3r::d3_nest(tm[, c("level1", "level2", "level3",
                      "level4", "level5", "level6",
                      "vSize")],
               value_cols = c("vSize"))


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
mapview::mapshot(p, url = file.path(plots_dir, "sample-distribution-pie.html"))
