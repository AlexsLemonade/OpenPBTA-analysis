# This script filters the given dataset to produce a summarized visualization
# of key variables within the dataset.
#
# Chante Bethell for CCDL 2019
#
# #### USAGE
# This script is intended to be run via the command line from the top directory
# of the repository as follows:
#
# Rscript analyses/sample-distribution-analysis/01-filter-across-types.R

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Function to filter based on primary_site
location_fn <- function(location) {
  # Given the name of a primary site, create a vector containing the disease
  # types expressed within the primary site.
  #
  # Note: the disease types found at all instances of the location substring
  #       will be included
  #
  # Args:
  #   location: the name of a primary site found within the dataset
  #
  # Returns:
  #   disease_type_vector: the vector of disease types expressed at the
  #                        named primary site
  disease_type_vector <- brain_location %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::filter(stringr::str_detect(primary_site, location)) %>%
    dplyr::pull(disease_type_new)
  unique(disease_type_vector)
}

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper execution, no matter where
# it is called from.
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Define the file.path to output directories
output_dir <- file.path(root_dir, "analyses", "sample-distribution-analysis")
results_dir <- file.path(output_dir, "results")
plots_dir <- file.path(output_dir, "plots")

# Create directories to hold the output.
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Read in dataset and remove NAs
histologies_df <-
  readr::read_tsv(file.path(root_dir, "data", "pbta-histologies.tsv")) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(disease_type_new))

# Filter the histologies file to account for multiple samples from the same
# individual and the fact that multiple experimental strategies are in this
# data.frame

# Retain only tumors for this analysis
histologies_df <- histologies_df %>%
  dplyr::filter(sample_type == "Tumor",
                composition == "Solid Tissue")

# data.frame with the count of each unique cancer type expression
disease_expression <- histologies_df %>%
  # some recurrences can have different disease_type_new values
  dplyr::distinct(Kids_First_Participant_ID, disease_type_new) %>%
  dplyr::group_by(disease_type_new) %>%
  dplyr::count(name = "count") %>%
  dplyr::arrange(dplyr::desc(count))

# Calculate the total count of the dataset
sum_count <- sum(disease_expression$count)

# Create a percent variable and round to 4 decimal places
# (so values will have 2 decimal places as percentages)
disease_expression <- disease_expression %>%
  dplyr::mutate(percent = paste0((round(count / sum_count, 4) * 100), "%"))

# Reorder the columns to be displayed in descending order by count on the plot
disease_expression$disease_type_new <- with(disease_expression,
                                            reorder(disease_type_new, -count))

# Write to tsv file
readr::write_tsv(disease_expression,
                 file.path(results_dir,
                           "disease_expression.tsv"))

# Create a bar plot of sample distribution across cancer types
gg_types <- disease_expression %>%
  ggplot2::ggplot(ggplot2::aes(x = disease_type_new, y = count, fill = count)) +
  ggplot2::geom_col() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Cancer Types", y = "Count",
                title = "Sample Distribution Across Cancer Types") +
  ggplot2::scale_y_continuous(breaks = seq(0, 500, by = 100)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 75,
    hjust = 1,
    size = 8
  ),
  panel.grid = ggplot2::element_blank()) +
  ggplot2::geom_text(nudge_y = 6.5, size = 2,
                     ggplot2::aes(label = paste0(disease_expression$percent)))

# Save plot
ggplot2::ggsave(
  gg_types,
  file = file.path(plots_dir, "distribution_across_cancer_types.pdf"),
  width = 22,
  height = 10
)

# data.frame with the location where each cancer type in the dataset is
# expressed, sorted to show highest expression
brain_location <- histologies_df %>%
  dplyr::distinct(Kids_First_Participant_ID, disease_type_new,
                  primary_site) %>%
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
                 file.path(results_dir, "primary_sites_counts.tsv"))
