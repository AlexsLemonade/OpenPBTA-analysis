# This script filters the given dataset to produce a summarized visualization
# of key variables within the dataset.
#
# Chante Bethell for CCDL 2019


# Load/install packages. This will be removed once a Dockerfile is established.
source("00-install-packages.R")

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Function to filter based on primary_site 
location_fn <- function(location) {
  disease_type_vector <- brain_location %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::filter(stringr::str_detect(primary_site, location)) %>%
    dplyr::pull(disease_type_new)
  unique(disease_type_vector)
}

# Create directories to hold the output.
if (!dir.exists("results")) {
  dir.create("results")
}
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Read in dataset
df2 <- data.frame(readr::read_tsv(
  file.path("..", "..", "data", "pbta-histologies.tsv")
))

# Remove na's
df2 <- df2 %>%
  dplyr::filter(!is.na(disease_type_new))

# data.frame with the count of each unique cancer type expression
disease_expression <- df2 %>%
  dplyr::group_by(disease_type_new) %>%
  dplyr::count(name = "count") %>%
  dplyr::arrange(dplyr::desc(count))

# Calculate the total count of the dataset
sum_count <- sum(disease_expression$count)

# Create a percent variable
disease_expression <- disease_expression %>%
  dplyr::mutate(percent = paste0(((count / sum_count) * 100), "%"))

# Reorder the columns to be displayed in descending order by count on the plot
disease_expression$disease_type_new <- with(disease_expression,
                                            reorder(disease_type_new, -count))

# Write to tsv file
readr::write_tsv(disease_expression,
                 file.path("results",
                           "disease_expression.tsv"))

# Create a bar plot of sample distribution across cancer types
gg_types <- disease_expression %>%
  ggplot2::ggplot(ggplot2::aes(x = disease_type_new, y = count, fill = count)) +
  ggplot2::geom_col() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Cancer Types", y = "Count",
                title = "Sample Distribution Across Cancer Types") +
  ggplot2::scale_y_continuous(breaks = seq(0, 500, by = 100)) +
  ggplot2::theme(axis.text.x = element_text(
    angle = 75,
    hjust = 1,
    size = 8
  ),
  panel.grid = element_blank()) +
  ggplot2::geom_text(nudge_y = 6.5, size = 2,
                     ggplot2::aes(label = paste0(disease_expression$percent)))

# Save plot
ggplot2::ggsave(
  gg_types,
  file = file.path("plots", "distribution_across_cancer_types.pdf"),
  width = 22,
  height = 10
)

# data.frame with the location where each cancer type in the dataset is 
# expressed, sorted to show highest expression
brain_location <- df2 %>%
  dplyr::select(disease_type_new, primary_site) %>%
  dplyr::group_by(disease_type_new, primary_site) %>%
  dplyr::tally() %>%
  dplyr::arrange(dplyr::desc(n))

# Make a vector of primary sites 
primary_sites_vector <- c(
  "Basal Ganglia",
  "Brain Stem- Midbrain",
  "Brain Stem-Medulla",
  "Cerebellum",
  "Frontal Lobe",
  "Parietal Lobe",
  "Temporal Lobe",
  "Occipital Lobe",
  "Spinal Cord",
  "Ventricles",
  "Thalamus",
  "Cranial Nerves NOS",
  "Hippocampus",
  "Meninges",
  "Optic Pathway",
  "Pineal Gland",
  "Pons",
  "Skull",
  "Spine NOS",
  "Suprasellar/Hypothalamic/Pituitary"
)

# This step helps us with melting
names(primary_sites_vector) <- primary_sites_vector

# For each string in primary sites vector use location_fn to get the vector of 
# disease types it's sorted by 
cancer_types_list <- lapply(primary_sites_vector, location_fn)

# Count the disease types for each primary site by taking the length of each 
# element of the list
cancer_types_counts <- lapply(cancer_types_list, length)

# Turn the list of counts into a data.frame and relabel, reorder columns
primary_sites_counts <- reshape2::melt(cancer_types_counts,
                                       value.name = "number_of_types") %>%
  dplyr::rename(primary_site = L1) %>%
  dplyr::select(primary_site, number_of_types) %>%
  # Get the max_type, it is the first element of the cancer_types_list.
  # Get the second_max_type, it is the second element.
  dplyr::mutate(max_type = unlist(lapply(cancer_types_list,
                                         function(x) x[1])),
                second_max_type = unlist(lapply(cancer_types_list,
                                                function(x) x[2]))) %>%
  # Reorder the rows
  dplyr::arrange(dplyr::desc(number_of_types))

# Write to tsv file
readr::write_tsv(primary_sites_counts,
                 file.path("results", "primary_sites_counts.tsv"))
sessionInfo()
