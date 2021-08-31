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
            "histology_label_color_table.tsv"),
  col_types = readr::cols())

# Read binary color for Male/Female color coding
binary_mapping <- readr::read_tsv(
  file.path(root_dir,
            "figures",
            "palettes",
            "binary_color_palette.tsv")
)

# Read tumor_descriptor palette
tumor_descriptor_palette <- readr::read_tsv(
  file.path(root_dir,
            "figures",
            "palettes",
            "tumor_descriptor_palette.tsv")
)

# Create final data.frame prepped for treemap and sunburst functions
final_df <- histologies_df %>%
  dplyr::filter(sample_type == "Tumor",
                composition == "Solid Tissue") %>%
  # Join on the color codes
  dplyr::inner_join(histology_label_mapping,
                    c("Kids_First_Biospecimen_ID", "sample_type",
                      "integrated_diagnosis", "Notes","harmonized_diagnosis",
                      "broad_histology", "short_histology", "cancer_group")) %>%
  # Extract WGS, WXS, RNA-Seq
  # To-do targeted sequencing?
  reshape2::dcast(sample_id + 
                    broad_histology +
                    cancer_group + 
                    tumor_descriptor + 
                    germline_sex_estimate +
                    cancer_group_hex_codes ~ experimental_strategy,
                  fun.aggregate = function(x){as.integer(length(x)>0)}) %>%
  dplyr::mutate(`DNA-Seq` = dplyr::if_else(WGS==1,"WGS","Not Available"),
                `DNA-Seq` = dplyr::if_else(WXS==1,"WXS",`DNA-Seq`),
                `DNA-Seq` = dplyr::if_else(`Targeted Sequencing`==1,"Targeted Sequencing",`DNA-Seq`),
                `RNA-Seq`= dplyr::if_else(`RNA-Seq`==1,"RNA-Seq","Not Available"),
                cancer_group = dplyr::if_else(is.na(cancer_group),
                                              "Other",cancer_group)) %>%
  # Get distinct based on participant IDs
  dplyr::distinct(sample_id, 
                  `DNA-Seq`,
                  `RNA-Seq`,
                  broad_histology,
                  cancer_group, 
                  tumor_descriptor,
                  germline_sex_estimate, 
                  cancer_group_hex_codes) %>% 
  # Select our columns of interest
  dplyr::select(broad_histology,
                cancer_group,
                `DNA-Seq`,
                `RNA-Seq`, 
                tumor_descriptor,
                germline_sex_estimate,
                cancer_group_hex_codes) %>%
  # Remove any row that has an NA
  dplyr::filter(complete.cases(.)) %>%
  # Group by all columns in order to count
  dplyr::group_by( broad_histology,
                   cancer_group,
                   `DNA-Seq`,
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
                level3 = `DNA-Seq`,
                level4 = `RNA-Seq`,
                level5 = tumor_descriptor,
                level6 = germline_sex_estimate) %>%
  # Reorder the rows according to the 3 levels
  dplyr::arrange(level1, level2, level3, level4, level5, level6) %>%
  dplyr::ungroup() %>%
  # tbl_df -> data.frame
  as.data.frame() 

# Save to tsv file
readr::write_tsv(final_df, file.path(results_dir, "plots_df.tsv"))

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
# merge to get hex_codes 
level1 <- tm %>%
  dplyr::filter(level == 1) %>%
  dplyr::select(-color)%>%
  dplyr::inner_join(histology_label_mapping
                          , by= c( "level1"="broad_histology")) %>%
  dplyr::rename(color=hex_codes)%>%
  dplyr::select(colnames(tm))

level2 <- tm %>%
  dplyr::filter(level == 2) %>%
  dplyr::select(-color)%>%
  dplyr::inner_join(histology_label_mapping
                   , by= c( "level2"="cancer_group")) %>%
  dplyr::rename(color=cancer_group_hex_codes)%>%
  dplyr::select(colnames(tm))

level3 <- tm %>%
  dplyr::filter(level == 3 ) %>%
  dplyr::mutate(color = dplyr::case_when(level3 == "WGS" ~ "#E5E5E5",
                                         level3 == "WXS" ~ "#B3B3B3",
                                         level3 == "Targeted Sequencing" ~ "6b6b6b",
                                         TRUE ~"#FFFFFF"))
level4 <- tm %>%
  dplyr::filter(level == 4 ) %>%
  dplyr::mutate( color = dplyr::if_else(level5 == "RNA-Seq", "#CCCCCC", "#FFFFFF"))

level5 <-  tm %>%
  dplyr::filter(level == 5 ) %>%
  dplyr::select(-color)%>%
  dplyr::inner_join(tumor_descriptor_palette,
                    by= c( "level5"="color_names")) %>%
  dplyr::rename(color=hex_codes)%>%
  dplyr::select(colnames(tm))

level6 <-  tm %>%
  dplyr::filter(level == 6 ) %>%
  dplyr::mutate(color = dplyr::case_when(
    level6 == "Male" ~ binary_mapping$hex_codes[1],
    level6 == "Female" ~ binary_mapping$hex_codes[2],
    # na_color
    TRUE ~ binary_mapping$hex_codes[3])
    )

new.tm <- dplyr::bind_rows(level1, level2, level3, level4, level5, level6) %>%
  unique()

# Convert the new.tm data.frame into a d3.js hierarchy object which is needed
# for the sund2b plot
tmnest <-
  d3r::d3_nest(new.tm,
               value_cols = colnames(new.tm)[-(1:6)])


# Create a sunburst plot
sun_plot <-
  sunburstR::sunburst(
    data = tmnest,
    valueField = "vSize",
    count = TRUE,
    sumNodes = FALSE
  )

# Create an interactive sund2b plot
p <- sunburstR::sund2b(tmnest,
                       colors = htmlwidgets::JS("function(name, d){return d.color || '#ccc';}"),
                       valueField = "vSize")

# Create HTML outputs for the interactive plots
mapview::mapshot(p, url = file.path(plots_dir, "histology-pie.html"))
