# Create color palettes for each kind of data in OpenPBTA-analysis
# C. Savonen for ALSF - CCDL 
#
# Usage: 
#  Anywhere a plot is being made, source these TSV file and use the color palette for 
#  each appropriate data type. 
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
  dplyr::mutate(short_histology = as.character(tidyr::replace_na(short_histology, "none"))) %>% 
  dplyr::arrange(short_histology)

# Get unique histology categories
histologies <- unique(metadata$short_histology)

# Create named list that can be used with dplyr::recode
# These hex codes were retrieved from http://phrogz.net/css/distinct-colors.html
# With the settings on default for 37 colors. 
histology_color_key <- c("#f23d3d", "#731d1d", "#b38686", "#cc5c33", "#331c0d", 
"#ffb380", "#b25f00", "#f2d6b6", "#736556", "#ffaa00", "#4c3d00", "#e2f200", 
"#919926", "#d6f2b6", "#304d26", "#00f241", "#009929", "#698c7c", "#39e6c3", 
"#005359", "#263233", "#00c2f2", "#40a6ff", "#406280", "#0044ff", "#00144d", 
"#acbbe6", "#7373e6", "#3d0099", "#c200f2", "#917399", "#731d6d", "#f279da", 
"#cc0052", "#994d6b", "#4d2636", "#ffbfd9")

# Bring along histologies names
names(histology_color_key) <- histologies

# NA for histology we want to match the standard NA color
histology_color_key[names(histology_color_key) == "none"] <- na_color

# Structure this as a data.frame and save to TSV
data.frame(color_names = names(histology_color_key), 
           hex_codes = unlist(histology_color_key)) %>% 
  readr::write_tsv(file.path(output_dir, "histology_color_palette.tsv"))

# Example Usage: 
# # Add a column that has the color for each sample
# histology_sample_df <- metadata %>% 
#  dplyr::select(Kids_First_Biospecimen_ID, short_histology) %>%
#  dplyr::mutate(sample_color = dplyr::recode(short_histology, 
#                                             !!!histology_color_key)) 
#
### 3) A gradient color scale for numeric data.
gradient_col_palette <- c("#f7f7f7",
                          "#f7fcf5",
                          "#e5f5e0",
                          "#c7e9c0",
                          "#a1d99b",
                          "#74c476",
                          "#41ab5d",
                          "#238b45",
                          "#006d2c",
                          "#00441b", 
                          na_color)  

gradient_col_names <- c("gradient_0", 
                        "gradient_1",
                        "gradient_2",
                        "gradient_3",
                        "gradient_4",
                        "gradient_5",
                        "gradient_6",
                        "gradient_7",
                        "gradient_8",
                        "gradient_9", 
                        "na_color") 
# Format as data.frame
gradient_df <- data.frame(color_names = gradient_col_names, 
                          hex_codes = gradient_col_palette) %>%
  readr::write_tsv(file.path(output_dir, "gradient_color_palette.tsv"))

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
                           "#053061", 
                           na_color)

divergent_color_names <- c("divergent_low_5",
                           "divergent_low_4",
                           "divergent_low_3",
                           "divergent_low_2",
                           "divergent_low_1",
                           "divergent_neutral",
                           "divergent_high_1",
                           "divergent_high_2",
                           "divergent_high_3",
                           "divergent_high_4",
                           "divergent_high_5", 
                           "na_color") 

# Format as data.frame
divergent_df <- data.frame(color_names = divergent_color_names, 
                           hex_codes = divergent_col_palette) %>%
  readr::write_tsv(file.path(output_dir, "divergent_color_palette.tsv"))
  
## Example usage: 
# divergent_col_val <- c(seq(from = min(divergent_variable), 
#                           to = 0,
#                           length = 5), 
#                       0, 
#                       seq(from = 0, 
#                           to = max(divergent_variable),
#                           length = 5)
#
# col_fun <- circlize::colorRamp2(divergent_col_val, 
#                                 divergent_col_palette)

### 5) A binary color key which are the most extreme colors in the divergent color scale. 
binary_col_palette <- c("#2166ac", 
                        "#b2182b", 
                         na_color) 

binary_color_names <- c("binary_1", 
                        "binary_2", 
                        "na_color") 

# Format as data.frame
binary_df <- data.frame(color_names = binary_color_names, 
                        hex_codes = binary_col_palette) %>%
  readr::write_tsv(file.path(output_dir, "binary_color_palette.tsv"))

#### Quick little function for writing HEX table in README so I can copy and paste it
#color_format <- function(color_df) {
#  #Description: Provide color data.frame, this outputs the text to render swatches 
#  #in a GitHub markdown table. 
#  colors <- gsub("#", "", color_df$hex_codes)
#  paste0("<br>", color_df$color_names, 
#         ":![", colors, "](https://placehold.it/150x40/", colors, "/FFFFFF?text=", colors, ")")
#}
