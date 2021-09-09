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

### 2) A gradient color scale for numeric data.
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

### 3) A divergent color scale for bidirectional numeric data.
divergent_col_palette <- c("#053061",
                           "#2166ac",
                           "#4393c3",
                           "#92c5de",
                           "#d1e5f0",
                           "#f7f7f7",
                           "#fddbc7",
                           "#f4a582",
                           "#d6604d",
                           "#b2182b",
                           "#67001f",
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

### 4) A binary color key which are the most extreme colors in the divergent color scale. 
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

### 6) An oncoprint color key for each CNV, SNV and Fusion category.
oncoprint_col_palette <- c("#35978f",
                           "#000000",
                           "#56B4E9",
                           "#FFBBFF",
                           "#F0E442",
                           "#191970",
                           "#545454",
                           "#CAE1FF",
                           "#FFE4E1",
                           "#CC79A7",
                           "#56B4E9",
                           "#7B68EE",
                           "#00F021",
                           "#313695",
                           "#abd9e9",
                           "#c51b7d",
                           "#0072B2",
                           "#D55E00",
                           "#FF0000",
                           "#CD96CD")

oncoprint_color_names <- c("Missense_Mutation",
                           "Nonsense_Mutation",
                           "Frame_Shift_Del",
                           "Frame_Shift_Ins",
                           "Splice_Site",
                           "Translation_Start_Site",
                           "Nonstop_Mutation",
                           "In_Frame_Del",
                           "In_Frame_Ins",
                           "Stop_Codon_Ins",
                           "Start_Codon_Del",
                           "Fusion",
                           "Multi_Hit",
                           "Hom_Deletion",
                           "Hem_Deletion",
                           "amplification",
                           "loss",
                           "gain",
                           "High_Level_Gain",
                           "Multi_Hit_Fusion")

# Format as data.frame
oncoprint_df <- data.frame(color_names = oncoprint_color_names,
                           hex_codes = oncoprint_col_palette) %>%
  readr::write_tsv(file.path(output_dir, "oncoprint_color_palette.tsv"))

#### Quick function for writing HEX code table in figures/README.md

 swatches_table <- function(color_df) {
  # For a given color data.frame with columns `hex_codes` and `color_names`,
  # this function outputs the text to render swatches in a GitHub markdown table.

  # These urls need the `#`` dropped
  colors <- gsub("#", "", color_df$hex_codes)

  # Paste the url together
  swatches <- paste0("<br>", color_df$color_names,
                     ":![",
                     colors,
                     "](https://placehold.it/150x40/",
                     colors,
                     "/FFFFFF?text=",
                     colors, ")")
  # Cat this out
  cat(swatches)
 }
 
 # Format as data.frame
 tumor_descriptor_palette <- data.frame(color_names = c("Initial CNS Tumor",
                                                        "Progressive",
                                                        "Progressive Disease Post-Mortem",
                                                        "Recurrence",
                                                        "Second Malignancy"),
                                        hex_codes = c("#709AE1FF",
                                                      "#075149FF",
                                                      "#075149FF",
                                                      "#FD8CC1FF",
                                                      "#FD7446FF")) %>%
   readr::write_tsv(file.path(output_dir, "tumor_descriptor_palette.tsv"))
 
#
# Usage:
# Put whichever color palette you are updating into function.
# Copy the output from this function and paste it into the appropriate section
# of the table in the figures/README.md
#
# swatches_table(color_df)
