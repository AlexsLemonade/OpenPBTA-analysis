# Stephanie J. Spielman for CCDL, 2022

# This script creates three S3 panels:
#  CN status heatmap
#  CNV breaks across chromothripsis regions box/jitter plot
#  SV breaks across chromothripsis regions box/jitter plot


## Load libraries ------------------------
library(tidyverse)
library(ComplexHeatmap)

# set seed for jitter plots
set.seed(1979)


## Define paths  -----------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses")
figure_dir <- file.path(root_dir, "figures")
palettes_dir <- file.path(figure_dir, "palettes")
supp_figures_dir <- file.path(figure_dir, "pdfs", "supp", "figs3", "panels")
if(!dir.exists(supp_figures_dir)){
  dir.create(supp_figures_dir, recursive = TRUE)
}

## Define output files
cn_heatmap_file <- file.path(supp_figures_dir, "cn_status_heatmap.pdf")
chromo_cnv_file <- file.path(supp_figures_dir, "chromothripsis_cnv.pdf")
chromo_sv_file <- file.path(supp_figures_dir, "chromothripsis_sv.pdf")

# Zenodo CSV output directory and file paths
zenodo_tables_dir <- file.path(root_dir, "tables", "zenodo-upload")
figS3c_csv <- file.path(zenodo_tables_dir, "figure-S3c-data.csv.gz") # 15 MB without gz
figS3d_csv <- file.path(zenodo_tables_dir, "figure-S3d-data.csv") 
figS3e_csv <- file.path(zenodo_tables_dir, "figure-S3e-data.csv") 


## Figure S3C ----------------------------
# Adapted from https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/cnv-chrom-plot/cn_status_heatmap.Rmd

# Define paths for this figure
original_dir <- file.path(analyses_dir, "cnv-chrom-plot")
cnv_dir <- file.path(analyses_dir, "copy_number_consensus_call", "results")

# Import bin calls
bin_calls_df <- read_tsv(file.path(original_dir, "results", "cn_status_bp_per_bin.tsv"))

# Import color palettes
histology_label_mapping <- readr::read_tsv(file.path(palettes_dir, "broad_histology_cancer_group_palette.tsv"))

# array of the three hex codes
divergent_col_hex <- readr::read_tsv(
  file.path(figure_dir, "palettes", "divergent_color_palette.tsv")
  ) %>%
  # Only keep only these three colors - we only need three values for `gain`, `neutral`, and `loss`.
  filter(
    color_names %in% c("divergent_low_4", "divergent_neutral", "divergent_high_4")
  ) %>%
  pull(hex_codes)


# Read in metadata and merge with colors
metadata <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"), guess_max = 10000) %>%
  # Join on the colors
  left_join(histology_label_mapping, by = c("broad_histology", "cancer_group")) %>%
  # Select what is needed
  select(contains("cancer_group"),
         Kids_First_Biospecimen_ID) %>%
  # remove "Other"
  filter(cancer_group != "Other")


### Set up chromosomal sizes for making bins.
chr_sizes <- read_tsv(
  file.path(data_dir, "WGS.hg38.strelka2.unpadded.bed"),
  col_names = c("chrom", "start", "end")
  ) %>%
  # Reformat the chromosome variable to drop the "chr"
  mutate(chrom = factor(gsub("chr", "", chrom),
                        levels = c(1:22, "X", "Y", "M"))
  ) %>%
  # Remove sex chromosomes
  filter(!(chrom %in% c("X", "Y", "M")))





########### Set up heatmap annotation objects


# Make chromosome size named vector for Heatmap annotation
chr_sizes_vector <- chr_sizes$end
names(chr_sizes_vector) <- chr_sizes$chrom

# Get a vector of the biospecimen IDs we use in the heatmap
sample_ids <- bin_calls_df$biospecimen_id


# Make color key.
color_key <- structure(
  c(divergent_col_hex, "#9932CC", "#fed8b1"),
  names = c("loss", "neutral", "gain", "unstable", "uncallable")
)


### Make chromosome labeling `HeatmapAnnotation` object

# First, set up bins of ~1Mb size and uncompress
bins <- GenomicRanges::tileGenome(
  chr_sizes_vector,
  tilewidth = 1e6
)
bins <- unlist(bins)
# use to set up chrs
chrs <- paste0("chr", S4Vectors::decode(bins@seqnames))
chrs <- factor(chrs, levels = paste0("chr", 1:22))

# Make a key for assigning alternating colors to the chromosomes
chr_colors <- rep(c("white", "black"), length.out = length(unique(chrs)))
names(chr_colors) <- unique(chrs)

# Get coordinate start, end, midpoint positions
chr_start <- match(unique(chrs), chrs)
chr_end <- chr_start + summary(chrs)
mid_points_chr <- floor((chr_start + chr_end) / 2)

# Make text labels for chromosome text
chr_text <- ComplexHeatmap::anno_mark(
  at = mid_points_chr,
  labels = levels(chrs),
  which = "column",
  side = "bottom",
  labels_rot = 45,
  labels_gp = grid::gpar(cex = 0.65)
)

# Create the Heatmap annotation object
chr_annot <- HeatmapAnnotation(
  df = data.frame(chrs),
  col = list(chrs = chr_colors),
  name = "",
  show_legend = FALSE,
  show_annotation_name = FALSE,
  mark = chr_text, # Put the text in
  border = TRUE
)



### Make row annotation object
# Get the histologies for the samples in this set and order them by histology

# Figure out cancer group counts, where NA display groups are _excluded_
cancer_group_counts <- metadata %>%
  drop_na(cancer_group_display) %>%
  # IDs of interest
  filter(Kids_First_Biospecimen_ID %in% sample_ids) %>%
  # count and create new label. we run factor code to ensure the right order
  count(cancer_group_display, name = "cancer_group_n") %>%
  mutate(
    cancer_group_display = fct_reorder(cancer_group_display, cancer_group_n, .desc=T),
    # and now make the label
    cancer_group_display_n = glue::glue("{cancer_group_display} (N = {cancer_group_n})")
  ) %>%
  # arrange so we're in the right order
  arrange(cancer_group_display)

samples_for_heatmap <- metadata %>%
  # again, filter to IDs of interest
  filter(Kids_First_Biospecimen_ID %in% sample_ids) %>%
  select(Kids_First_Biospecimen_ID,
         cancer_group_display, # display names needed for factoring
         cancer_group_hex) %>%
  inner_join(
    select(
      cancer_group_counts,
      cancer_group_n,
      cancer_group_display
    )
  ) %>%
  mutate(
    cancer_group_display = fct_reorder(cancer_group_display, cancer_group_n, .desc=T),
    # and now relabel
    cancer_group_display = factor(cancer_group_display, labels = cancer_group_counts$cancer_group_display_n)
  ) %>%
  # Make sure the rows are also in cancer_group_order
  arrange(cancer_group_display) %>%
  # Store as rownames
  column_to_rownames("Kids_First_Biospecimen_ID") %>%
  # remove coumn
  select(-cancer_group_n)

# Make a color key that's formatted for ComplexHeatmap
# Get a distinct version of the color keys
cancer_color_key_df <- samples_for_heatmap %>%
  dplyr::select(cancer_group_display, cancer_group_hex) %>%
  dplyr::distinct()

# Make color key specific to these samples
cancer_color_key <- cancer_color_key_df$cancer_group_hex
names(cancer_color_key) <- cancer_color_key_df$cancer_group_display

# Get coordinate start positions
hist_start <- match(names(cancer_color_key), samples_for_heatmap$cancer_group_display)

# Get coordinate end positions for each histology group
hist_end <- hist_start + summary(samples_for_heatmap$cancer_group_display)

# Get mid points of histology group for labels
mid_points_hist <- floor((hist_start + hist_end) /2)


# Make text labels for chromosome text
hist_text <- ComplexHeatmap::anno_mark(
  at = mid_points_hist,
  labels = levels(samples_for_heatmap$cancer_group_display),
  which = "row",
  side = "right",
  labels_gp = grid::gpar(cex = 0.5),
  link_width = grid::unit(5, "mm")
)

# Create the Heatmap annotation object
hist_annot <- ComplexHeatmap::HeatmapAnnotation(
  df = data.frame(samples_for_heatmap %>% select(-cancer_group_hex)), # select out to avoid a double!
  col = list(cancer_group_display = cancer_color_key),
  which = "row",
  show_annotation_name = FALSE,
  show_legend = FALSE,
  mark = hist_text, # Put the text in
  border = TRUE
  )


# Format `bin_calls_df` as a matrix with rownames for `ComplexHeatmap` to use.
bin_calls_mat <- bin_calls_df %>%
  tibble::column_to_rownames("biospecimen_id") %>%
  as.matrix()
# Ensure that this matrix is in the same order as the annotation and double check order
bin_calls_mat <- bin_calls_mat[rownames(samples_for_heatmap), ]
if (all.equal(rownames(bin_calls_mat), rownames(samples_for_heatmap)) != TRUE) {
  stop("Bad data processing for CN status heatmap")
}

## Assemble CN status heatmap
heatmap <- ComplexHeatmap::Heatmap(
  bin_calls_mat,
  name = "CN status",
  col = color_key,
  row_split = samples_for_heatmap$cancer_group_display,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  rect_gp = grid::gpar(col = "black", lwd = .0000005),
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_labels = FALSE,
  bottom_annotation = chr_annot,
  right_annotation = hist_annot,
  heatmap_legend_param = list(nrow = 1, border = "black"),
  raster_quality = 8,
  border = "black",
  row_title = NULL
)

# Save heatmap.
pdf(cn_heatmap_file, width = 8, height = 7)
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

## Figure S3D-E ---------------------------
# originally from https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/chromothripsis/04-plot-chromothripsis-and-breakpoint-data.Rmd

# Set theme
theme_set(ggpubr::theme_pubr() + theme(legend.position = "none"))

# read in and prepare data
chromo_dir <- file.path(analyses_dir, "chromothripsis")
breakpoint_dir <- file.path(analyses_dir, "chromosomal-instability", "breakpoint-data")

chromoth_per_sample <- read_tsv(
  file.path(chromo_dir, "results", "chromothripsis_summary_per_sample.txt")
  )
cnv_densities <- read_tsv(
  file.path(breakpoint_dir, "cnv_breaks_densities.tsv")
  )
sv_densities <- read_tsv(
  file.path(breakpoint_dir, "sv_breaks_densities.tsv")
  )


# Rename columns before merging
cnv_densities <- cnv_densities %>%
  rename(Kids_First_Biospecimen_ID = samples,
         cnv_breaks_count = breaks_count)

sv_densities <- sv_densities %>%
  rename(Kids_First_Biospecimen_ID = samples) %>%
  rename(sv_breaks_count = breaks_count)

# Merge chromothripsis data and breakpoint data
merged_data <- chromoth_per_sample %>%
  inner_join(cnv_densities, by = "Kids_First_Biospecimen_ID") %>%
  inner_join(sv_densities, by = "Kids_First_Biospecimen_ID") %>%
  # Truncate # chromothripsis regions above 5
  mutate(count_regions_any_conf_truncated = ifelse(count_regions_any_conf>=5, ">=5", count_regions_any_conf),
         count_regions_any_conf_truncated = forcats::fct_relevel(count_regions_any_conf_truncated, ">=5", after = Inf))


##### Figure S3D
fig_s3d <- ggplot(merged_data) + 
  aes(x = count_regions_any_conf_truncated,
      y = cnv_breaks_count) +
  geom_jitter(width = 0.3, alpha = 0.5, size = 0.4) +
  geom_boxplot(color = "black", alpha = 0, outlier.shape=NA,
               # fatten controls *median line* width
               size = 0.4, fatten = 1) +
  xlab("Number of Chromothripsis Regions") +
  ylab("Number of CNV Breaks") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25)
  )

ggsave(chromo_cnv_file, fig_s3d, width = 2.65, height = 1.65,
       useDingbats=FALSE)

#### Figure S3E
fig_s3e <- ggplot(merged_data) + 
  aes(x = count_regions_any_conf_truncated,
      y = sv_breaks_count) +
  geom_jitter(width = 0.3, alpha = 0.5, size = 0.4) +
  geom_boxplot(color = "black", alpha = 0, outlier.shape=NA,
               # fatten controls *median line* width
               size = 0.4, fatten = 1) +
  xlab("Number of Chromothripsis Regions") +
  ylab("Number of SV Breaks") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25)
  )
ggsave(chromo_sv_file, fig_s3e, width = 2.65, height = 1.65,
       useDingbats=FALSE)


# Export CSVs for Zenodo upload ----------------

# Panel S3C

# To make output more informative, we'll change column names from just a bin number to a format:
#  `chrZ_binX`, where Z is the chromosome number and X is the bin number (the current column name)
bin_nums <- as.numeric(colnames(bin_calls_df)[-1])

# Create bin position names
new_bin_colnames <- glue::glue("{chrs}:{bins@ranges@start}-{bins@ranges@start + bins@ranges@width - 1}")
  
# Re-assign the bin column names, but don't change the sample column
colnames(bin_calls_df)[-1] <- new_bin_colnames
  

# Finally, we can export:
bin_calls_df %>%
  # rename to standardized name
  dplyr::rename(Kids_First_Biospecimen_ID = biospecimen_id) %>%
  # filter to samples that are actually in the plot
  dplyr::filter(Kids_First_Biospecimen_ID %in% rownames(samples_for_heatmap)) %>%
  # arrange on sample; note this column is already first
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  # export
  readr::write_csv(figS3c_csv)



# Panel S3D and S3E are made from the same data, 
#  so prep first and then subset to relevant variables
merged_data_export <- merged_data %>%
  # select just the columns in the plots
  dplyr::select(Kids_First_Biospecimen_ID, 
                count_regions_any_conf_truncated, 
                cnv_breaks_count, 
                sv_breaks_count) %>%
  # arrange on sample
  dplyr::arrange(Kids_First_Biospecimen_ID) 

# S3D
merged_data_export %>%
  # remove column that is in S3E
  select(-sv_breaks_count) %>%
  readr::write_csv(figS3d_csv)

# S3E
merged_data_export %>%
  # remove column that is in S3D
  select(-cnv_breaks_count) %>%
  readr::write_csv(figS3e_csv)

