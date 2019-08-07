# This script filters the given dataset to produce a summarized visualization of key variables within the dataset.
#
# CCDL - 2019 
# Chante Bethell 


# Load/install packages. This will be removed once a Dockerfile is established. 
source(file.path("..", "scripts", "install-packages.R"))

# magrittr pipe
`%>%` <- dplyr::`%>%`

# Function to filter dataset based on primary_site of disease, while arranging the cancer types in order of descending expression
location_fn <- function(location) {
  brain.location %>%
    arrange(desc(n)) %>%
    filter(str_detect(primary_site, location))
}

# Function to pad x with NA's
na.pad <- function(x, len){
  x[1:len]
}

# Function to make a padded data.frame
makePaddedDataFrame <- function(l, ...){
  maxlen <- max(sapply(l, length))
  data.frame(lapply(l, na.pad, len = maxlen), ...)
}

# Create directories to hold the output.
if (!dir.exists("../results")) {
  dir.create("../results")
}
if (!dir.exists("../plots")) {
  dir.create("../plots")
}

# Read in dataset 
df2 <- data.frame(read_tsv(file.path("..", "data", "2019-07-26-pbta-histologies - 2019-07-26-including-germline.tsv")))

# Remove na's 
df2 <- df2 %>%
  filter(!is.na(disease_type_new))

# data.frame with the count of each unique cancer type expression 
disease.expression <- df2 %>%
  group_by(disease_type_new) %>%
  filter(!is.na(disease_type_new)) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(count = n)

# Calculate the total count of the dataset 
sum.count <- sum(disease.expression$count)

# Create a percent variable 
disease.expression <- disease.expression %>%
  mutate(percent = (count/sum.count))

disease.expression$percent <- percent(disease.expression$percent)
disease.expression$disease_type_new <- with(disease.expression, reorder(disease_type_new, -count))

# Save the disease.expression data.frame as a table 
disease.expression.table <- kable(disease.expression, booktabs = TRUE) %>%
  kable_styling(bootstrap_options = "condensed", full_width = FALSE, 
                font_size = 6, latex_options = "scale_down") %>%
  save_kable(file.path("..", "results", "disease_expression_table.pdf"))

# Write to tsv file
write_tsv(disease.expression, file.path("..", "results", "disease_expression.csv"))

# Create a bar plot of sample distribution across cancer types 
gg.types <- disease.expression %>%
  ggplot(aes(x = disease_type_new, y = count, fill = count)) + 
  geom_col() +
  theme_bw() +
  labs(x = "Cancer Types", y = "Count", 
       title = "Sample Distribution Across Cancer Types") +
  scale_y_continuous(breaks = seq(0, 500, by = 100)) + 
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 8), 
        panel.grid = element_blank()) +
  geom_text(nudge_y = 6.5, size = 2, 
            aes(label = paste0(disease.expression$percent)))

# Save plot 
ggsave(gg.types, file = file.path("..", "plots", "distribution_across_cancer_types.pdf"), 
       width = 22, height = 10)

# data.frame with the location where each cancer type in the dataset is expressed sorted to show highest expression
brain.location <- df2 %>%
  select(disease_type_new, primary_site) %>%
  group_by(primary_site) %>%
  arrange(primary_site) %>%
  group_by(disease_type_new, primary_site) %>%
  tally() %>%
  arrange(desc(n)) %>%
  filter(!is.na(disease_type_new))


# Use the location_fn to create a data.frame for each unique primary site 
basal.ganglia <- location_fn("Basal Ganglia") 
basal.ganglia <- unique(basal.ganglia$disease_type_new)

brain.stem <- location_fn("Brain Stem")
brain.stem <- unique(brain.stem$disease_type_new)

cerebellum <- location_fn("Cerebellum")
cerebellum <- unique(cerebellum$disease_type_new)

frontal.lobe <- location_fn("Frontal Lobe")
frontal.lobe <- unique(frontal.lobe$disease_type_new)

parietal.lobe <- location_fn("Parietal Lobe")
parietal.lobe <- unique(parietal.lobe$disease_type_new)

temporal.lobe <- location_fn("Temporal Lobe")
temporal.lobe <- unique(temporal.lobe$disease_type_new)
  
occipital.lobe <- location_fn("Occipital Lobe")
occipital.lobe <- unique(occipital.lobe$disease_type_new)

spinal.cord <- location_fn("Spinal Cord")
spinal.cord <- unique(spinal.cord$disease_type_new)

ventricles <- location_fn("Ventricles")
ventricles <- unique(ventricles$disease_type_new)

thalamus <- location_fn("Thalamus")
thalamus <- unique(thalamus$disease_type_new)

cranial.nerves <- location_fn("Cranial Nerves NOS")
cranial.nerves <- unique(cranial.nerves$disease_type_new)

hippocampus <- location_fn("Hippocampus")
hippocampus <- unique(hippocampus$disease_type_new)

meninges <- location_fn("Meninges")
meninges <- unique(meninges$disease_type_new)

optic.pathway <- location_fn("Optic Pathway")
optic.pathway <- unique(optic.pathway$disease_type_new)

pineal.gland <- location_fn("Pineal Gland")
pineal.gland <- unique(pineal.gland$disease_type_new)

suprasellar.hypothalamic.pituitary <- location_fn("Suprasellar/Hypothalamic/Pituitary")
suprasellar.hypothalamic.pituitary <- unique(suprasellar.hypothalamic.pituitary$disease_type_new)

pons <- location_fn("Pons")
pons <- unique(pons$disease_type_new)

skull <- location_fn("Skull")
skull <- unique(skull$disease_type_new)

spine <- location_fn("Spine NOS")
spine <- unique(spine$disease_type_new)

# Make a data.frame containing each unique primary site and the cancer types therein 
primary.sites <- makePaddedDataFrame(list(frontal.lobe = frontal.lobe, parietal.lobe = parietal.lobe, 
                                          temporal.lobe = temporal.lobe, spinal.cord = spinal.cord,
                                          ventricles = ventricles, cerebellum = cerebellum, 
                                          occipital.lobe = occipital.lobe, 
                                          suprasellar.hypothalamic.pituitary = suprasellar.hypothalamic.pituitary,
                                          skull = skull, brain.stem = brain.stem, thalamus = thalamus, 
                                          spine = spine, pons = pons,
                                          pineal.gland = pineal.gland, meninges = meninges,
                                          basal.ganglia = basal.ganglia, optic.pathway = optic.pathway,
                                          hippocampus = hippocampus, cranial.nerves = cranial.nerves
                                          ))
# Save table to pdf 
primary.site.table <- kable(primary.sites, booktabs = TRUE) %>%
  kable_styling(bootstrap_options = "condensed", full_width = FALSE, font_size = 6, latex_options = "scale_down") %>%
  save_kable(file.path("..", "results", "primary_sites_table.pdf"))

# Write to tsv file
write_tsv(primary.sites, file.path("..", "results", "primary_sites.tsv"))

sessionInfo()
