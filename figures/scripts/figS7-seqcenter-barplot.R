# S Spielman for CCDL, 2023
# Create panel A for Figure S7 that shows how histologies (cancer groups)
# are not balanced across RNA library preparations (polyA vs. stranded)


#### Libraries -----------------------------------------------------------------

library(tidyverse)
library(ggpubr)

#### Directories and files -----------------------------------------------------

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
figures_dir <- file.path(root_dir, "figures")
output_dir <- file.path(figures_dir, "pdfs", "supp", "figS7", "panels")

output_pdf <- file.path(output_dir, "RNA_library_barplot.pdf")

# Zenodo CSV output directory and file path
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
figS7a_csv <- file.path(zenodo_tables_dir, "figure-S7a-data.csv")


metadata_file <- file.path(data_dir, "pbta-histologies.tsv")
palette_file <- file.path(figures_dir, "palettes", "broad_histology_cancer_group_palette.tsv")
binary_pal_file <- file.path(figures_dir, "palettes", "binary_color_palette.tsv")

##### Read in and prepare data for plot ----------------------

metadata_df <- read_tsv(metadata_file, guess_max = 10000)
palette_df <- read_tsv(palette_file)
binary_pal_df <- read_tsv(binary_pal_file)

# Get the RNA samples and combine with `cancer_group_display` information
library_cancer_group_df <- metadata_df %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  # We want one row per patient here, using `sample_id`
  mutate(dup_participant = duplicated(sample_id)) %>%
  filter(!dup_participant) %>%
  # Keep RNA_library and cancer group for each
  select(Kids_First_Biospecimen_ID, RNA_library, cancer_group) %>%
  # Cancer group display
  inner_join(
    select(palette_df,
           cancer_group,
           cancer_group_display),
    by = "cancer_group"
  ) %>%
  # Remove NAs
  drop_na(cancer_group_display) %>%
  # Arrange cancer_group_display factor in order of total counts,
  #  making sure "Other" is still last
  mutate(cancer_group_display = stringr::str_wrap(cancer_group_display, 25), # make the labels fit
         cancer_group_display = forcats::fct_infreq(cancer_group_display),
         cancer_group_display = forcats::fct_relevel(cancer_group_display,
                                                     "Other", after = Inf)
  )

# Prepare binary color palette
pal_colors <- rev(binary_pal_df$hex_codes[binary_pal_df$color_names != "na_color"])

###### Make the plot --------------------------

rna_library_barplot <- ggplot(library_cancer_group_df) +
  aes(x = cancer_group_display,
      fill = RNA_library) +
  geom_bar(color = "black", size = 0.25) +
  scale_fill_manual(name = "RNA library preparation",
                    values = pal_colors) +
  labs(
    x = "Cancer group",
    y = "Total number of participants with RNA-Seq"
  ) +
  ggpubr::theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10)
  )

ggsave(
  output_pdf,
  rna_library_barplot,
  width = 8.5,
  height = 6
)


# Export CSV for Zenodo upload
# no samples so nothing to arrange
library_cancer_group_df %>%
  # remove \n
  dplyr::mutate(cancer_group_display = stringr::str_replace(cancer_group_display, "\n", " ")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cancer_group_display, RNA_library) %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>% 
  readr::write_csv(figS7a_csv)


