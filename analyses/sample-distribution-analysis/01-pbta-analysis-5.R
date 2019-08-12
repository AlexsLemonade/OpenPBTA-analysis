# This script filters the given dataset to produce a summarized visualization
# of key variables within the dataset.
#
# Chante Bethell for CCDL 2019


# Load/install packages. This will be removed once a Dockerfile is established.
source("00-install-packages.R")

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Function to filter dataset based on primary_site of disease, while arranging
# the cancer types in order of descending expression.
location_fn <- function(location) {
  brain_location %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::filter(stringr::str_detect(primary_site, location))
}

# Function to pad x with NA's
na_pad <- function(x, len) {
  x[1:len]
}

# Function to make a padded data.frame
make_padded_data_frame <- function(l, ...) {
  maxlen <- max(sapply(l, length))
  data.frame(lapply(l, na_pad, len = maxlen), ...)
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
  file.path("..", "..", "data",
            "2019-07-26-pbta-histologies - 2019-07-26-including-germline.tsv"
  )
))

# Remove na's
df2 <- df2 %>%
  dplyr::filter(!is.na(disease_type_new))

# data.frame with the count of each unique cancer type expression
disease_expression <- df2 %>%
  dplyr::group_by(disease_type_new) %>%
  dplyr::filter(!is.na(disease_type_new)) %>%
  dplyr::count() %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::rename(count = n)

# Calculate the total count of the dataset
sum_count <- sum(disease_expression$count)

# Create a percent variable
disease_expression <- disease_expression %>%
  dplyr::mutate(percent = (count / sum_count))

# Format the percent values to include %
disease_expression$percent <-
  formattable::percent(disease_expression$percent)

# Reorder the columns to be displayed in descending order by count on the plot
disease_expression$disease_type_new <- with(disease_expression,
                                            reorder(disease_type_new, -count))

# Write to tsv file
readr::write_tsv(disease_expression,
                 file.path("results",
                           "disease_expression.tsv"))

# Create a bar plot of sample distribution across cancer types
gg_types <- disease_expression %>%
  ggplot2::ggplot(aes(x = disease_type_new, y = count, fill = count)) +
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
                     aes(label = paste0(disease_expression$percent)))

# Save plot
ggplot2::ggsave(
  gg_types,
  file = file.path("plots", "distribution_across_cancer_types.pdf"),
  width = 22,
  height = 10
)

# data.frame with the location where each cancer type in the dataset is expressed 
# sorted to show highest expression
brain_location <- df2 %>%
  dplyr::select(disease_type_new, primary_site) %>%
  dplyr::group_by(primary_site) %>%
  dplyr::arrange(primary_site) %>%
  dplyr::group_by(disease_type_new, primary_site) %>%
  dplyr::tally() %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(!is.na(disease_type_new))


# Use the location_fn to create a data.frame for each unique primary site
basal_ganglia <- location_fn("Basal Ganglia")
basal_ganglia <- unique(basal_ganglia$disease_type_new)

brain_stem <- location_fn("Brain Stem")
brain_stem <- unique(brain_stem$disease_type_new)

cerebellum <- location_fn("Cerebellum")
cerebellum <- unique(cerebellum$disease_type_new)

frontal_lobe <- location_fn("Frontal Lobe")
frontal_lobe <- unique(frontal_lobe$disease_type_new)

parietal_lobe <- location_fn("Parietal Lobe")
parietal_lobe <- unique(parietal_lobe$disease_type_new)

temporal_lobe <- location_fn("Temporal Lobe")
temporal_lobe <- unique(temporal_lobe$disease_type_new)

occipital_lobe <- location_fn("Occipital Lobe")
occipital_lobe <- unique(occipital_lobe$disease_type_new)

spinal_cord <- location_fn("Spinal Cord")
spinal_cord <- unique(spinal_cord$disease_type_new)

ventricles <- location_fn("Ventricles")
ventricles <- unique(ventricles$disease_type_new)

thalamus <- location_fn("Thalamus")
thalamus <- unique(thalamus$disease_type_new)

cranial_nerves <- location_fn("Cranial Nerves NOS")
cranial_nerves <- unique(cranial_nerves$disease_type_new)

hippocampus <- location_fn("Hippocampus")
hippocampus <- unique(hippocampus$disease_type_new)

meninges <- location_fn("Meninges")
meninges <- unique(meninges$disease_type_new)

optic_pathway <- location_fn("Optic Pathway")
optic_pathway <- unique(optic_pathway$disease_type_new)

pineal_gland <- location_fn("Pineal Gland")
pineal_gland <- unique(pineal_gland$disease_type_new)

pons <- location_fn("Pons")
pons <- unique(pons$disease_type_new)

skull <- location_fn("Skull")
skull <- unique(skull$disease_type_new)

spine <- location_fn("Spine NOS")
spine <- unique(spine$disease_type_new)

suprasellar_hypothalamic_pituitary <-
  location_fn("Suprasellar/Hypothalamic/Pituitary")
suprasellar_hypothalamic_pituitary <-
  unique(suprasellar_hypothalamic_pituitary$disease_type_new)

# data.frame containing each unique primary site and the cancer types therein
primary_sites <-
  make_padded_data_frame(
    list(
      frontal_lobe = frontal_lobe,
      parietal_lobe = parietal_lobe,
      temporal_lobe = temporal_lobe,
      spinal_cord = spinal_cord,
      ventricles = ventricles,
      cerebellum = cerebellum,
      occipital_lobe = occipital_lobe,
      suprasellar_hypothalamic_pituitary = suprasellar_hypothalamic_pituitary,
      skull = skull,
      brain_stem = brain_stem,
      thalamus = thalamus,
      spine = spine,
      pons = pons,
      pineal_gland = pineal_gland,
      meninges = meninges,
      basal_ganglia = basal_ganglia,
      optic_pathway = optic_pathway,
      hippocampus = hippocampus,
      cranial_nerves = cranial_nerves
    )
  )

# Write to tsv file
readr::write_tsv(primary_sites, file.path("results", "primary_sites.tsv"))

# data.frame to put together the total number of types per site
primary_sites_counts <-
  data.frame(apply(primary_sites, 2, function(x)
    length(which(!is.na(x))))) %>%
  dplyr::rename(number_of_types = apply.primary_sites..2..function.x..length.which..is.na.x....) %>%
  tibble::rownames_to_column("primary_site")

# Add max_type column to give the highest expressed type at each site
primary_sites_counts$max_type <- as.factor(t(primary_sites[1, ]))

# Add second_max_type column to give the second highest expressed type at 
# each site. This was done to represent a broader scope of types. 
primary_sites_counts$second_max_type <- as.factor(t(primary_sites[2,]))

# Write to tsv file
readr::write_tsv(primary_sites_counts,
                 file.path("results", "primary_sites_counts.tsv"))

sessionInfo()
