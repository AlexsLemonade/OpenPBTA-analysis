# Create color palettes for each kind of data in OpenPBTA-analysis
# C. Savonen for ALSF - CCDL 
#
# Usage: 
#  Anywhere a plot is being made, source this file and use the color palette for 
#  each appropriate data type. 
#
# Example: 
# source(file.path("scripts", "color_palettes.R"))
#
#
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Output to palette directory
output_dir <- file.path(root_dir, "figures", "palettes")

### 1) A color for NA data
na_color <- "#f1f1f1"

### 2) Color palette for each histology group in short_histology.

# Read in metadata file
metadata <- readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv")) %>% 
  dplyr::mutate(short_histology = as.character(tidyr::replace_na(short_histology, "none")))

# Get unique histology categories
histologies <- unique(metadata$short_histology)

# Get number of unique histology categories
n_histologies <- length(histologies)

# Create named list that can be used with dplyr::recode
histologies_color_key <- hsv(h = 1:n_histologies/n_histologies * .85, 
                             v = c(.8,1,1), 
                             s = c(1,1, .6))

# Bring along histologies names
names(histologies_color_key) <- histologies

# NA for histology we want to match the standard NA color
histologies_color_key[names(histologies_color_key) == "none"] <- na_color

# Add a column that has the color for each sample
histology_col_by_sample <- metadata %>% 
  dplyr::select(Kids_First_Biospecimen_ID, short_histology) %>%
  dplyr::mutate(sample_color = dplyr::recode(short_histology, 
                                             !!!histologies_color_key)) %>% 
  readr::write_tsv(file.path(output_dir, "histology_color_by_sample.tsv"))

# Drop "none" category from histology color key
histologies_color_key <- histologies_color_key[histologies_color_key != "none"]

# Example Usage: 
# histologies_w_color_key <-
#  data.frame(Kids_First_Biospecimen_ID = common_samples) %>%
#  dplyr::inner_join(histology_col_by_sample)
#
### 3) A gradient color scale for numeric data.
gradient_col_palette <- c("#f7fcf5",
                          "#e5f5e0",
                          "#c7e9c0",
                          "#a1d99b",
                          "#74c476",
                          "#41ab5d",
                          "#238b45",
                          "#006d2c",
                          "#00441b")

## Example usage for variable: 
# gradient_col_val <- seq(from = min(variable), to = max(variable), 
#               length.out = length(gradient_col_palette))
#
# col_fun <- circlize::colorRamp2(gradient_col_val, 
#                                 gradient_col_palette)

### 4) A divergent color scale for bidirectional numeric data. 
divergent_col_palette <- c("#67001f",
                           "#b2182b",
                           "#d6604d",
                           "#f4a582",
                           "#fddbc7",
                           "#f7f7f7",
                           "#d1e5f0",
                           "#92c5de",
                           "#4393c3",
                           "#2166ac",
                           "#053061")
## Example usage: 
# gradient_col_val <- c(seq(from = min(divergent_variable), 
#                           to = 0,
#                           length = 5), 
#                       0, 
#                       seq(from = 0, 
#                           to = max(divergent_variable),
#                           length = 5)
#
# col_fun <- circlize::colorRamp2(gradient_col_val, 
#                                 gradient_col_palette)

### 5) A binary color key which are the most extreme colors in the divergent color scale. 
binary_col_palette <- c("#67001f", 
                        "#053061")
## Example usage: 
# names(binary_col_palette) <- c("Amp", "Del")

######################## Write these to an RDS file ############################
readr::write_rds(list(na_color,
                      histologies_color_key,
                      gradient_col_palette, 
                      divergent_col_palette, 
                      binary_col_palette),
                 file.path(output_dir, "hex_color_palettes.rds"))
