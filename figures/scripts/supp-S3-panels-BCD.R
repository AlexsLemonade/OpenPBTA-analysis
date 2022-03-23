# Stephanie J. Spielman for CCDL, 2022

# This script creates three S3 panels (B, C, D)


## Load libraries ------------------------
library(tidyverse)
library(ComplexHeatmap)


## Define paths  -----------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses")
scratch_dir <- file.path(root_dir, "scratch")
figure_dir <- file.path(root_dir, "figures")
palettes_dir <- file.path(figure_dir, "palettes")
supp_figures_dir <- file.path(figure_dir, "pdfs", "supp")

## Define output files
figure_s3b_file <- file.path(supp_figures_dir, "supp_figS3-B.pdf")
figure_s3c_file <- file.path(supp_figures_dir, "supp_figS3-C.pdf")
figure_s3d_file <- file.path(supp_figures_dir, "supp_figS3-D.pdf")

## Figure S3B ----------------------------
# Adapted from https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/cnv-chrom-plot/cn_status_heatmap.Rmd

# Source custom functions and define paths for this figure
original_dir <- file.path(analyses_dir, "cnv-chrom-plot")
cnv_dir <- file.path(analyses_dir, "copy_number_consensus_call", "results")

source(file.path(original_dir, "util", "bin-coverage.R"))

# Set cutoffs
length_max <- 1e7 # The max length of a segment to use the data.
frac_uncallable <- 0.75 # Set minimum percentage of a bin that should be callable to report data.
frac_threshold <- 0.75 # Absolute fraction needed for a bin to be called a particular status
min_group_size <- 2 # Any groups smaller than this will be added into the `Other` group for the resulting heatmap

# Import color palettes
histology_label_mapping <- readr::read_tsv(file.path(palettes_dir, "broad_histology_cancer_group_palette.tsv"))# %>% 
 # select(contains("cancer_group")) 

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
  left_join(histology_label_mapping, by = "cancer_group") %>%
  # Select what is needed
  select(contains("cancer_group"), 
                Kids_First_Biospecimen_ID, 
                tumor_ploidy)


### Set up consensus copy number data
# Read in the segment copy number data
seg_data <- data.table::fread(
  file.path(
    cnv_dir,
    "pbta-cnv-consensus.seg.gz"
  ),
  data.table = FALSE
)

# Set up the status for each consensus segment. 
seg_data <- seg_data %>%
  # Join the histology column to this data
  inner_join(
    select(
      metadata,
      Kids_First_Biospecimen_ID,
      cancer_group_display,
      tumor_ploidy
    ),
    by = c("ID" = "Kids_First_Biospecimen_ID")
  ) %>%
  # Reformat the chromosome variable to drop the "chr"
  mutate(chrom = factor(gsub("chr", "", chrom),
                               levels = c(1:22, "X", "Y")
  )) %>%
  # Recode the copy number status based on ploidy
  mutate(status = case_when(
    # when the copy number is less than inferred ploidy, mark this as a loss
    copy.num < tumor_ploidy ~ "loss",
    # if copy number is higher than ploidy, mark as a gain
    copy.num > tumor_ploidy ~ "gain",
    copy.num == tumor_ploidy ~ "neutral"
  )) %>%
  # Remove sex chromosomes
  filter(
    !(chrom %in% c("X", "Y", "M")),
    !is.na(status)
  )


# Set up seg data as GenomicRanges. 
seg_ranges <- GenomicRanges::GRanges(
  seqnames = seg_data$chrom,
  ranges = IRanges::IRanges(
    start = seg_data$loc.start,
    end = seg_data$loc.end
  ),
  status = seg_data$status,
  histology = seg_data$cancer_group,
  biospecimen = seg_data$ID
)

# Filter out segments that are longer than our cutoff. 
filtered_seg_ranges <- seg_ranges[which(seg_ranges@ranges@width < length_max)]


### Set up chromosomal sizes for making bins. 
chr_sizes <- read_tsv(file.path(data_dir, "WGS.hg38.strelka2.unpadded.bed"),
                             col_names = c("chrom", "start", "end")
) %>%
  # Reformat the chromosome variable to drop the "chr"
  mutate(chrom = factor(gsub("chr", "", chrom),
                               levels = c(1:22, "X", "Y", "M")
  )) %>%
  # Remove sex chromosomes
  filter(!(chrom %in% c("X", "Y", "M")))


# Make chromosome size named vector for Heatmap annotation
chr_sizes_vector <- chr_sizes$end
names(chr_sizes_vector) <- chr_sizes$chrom


### Set up uncallable regions data 
uncallable_bed <- readr::read_tsv(
  file.path(
    analyses_dir,
    "copy_number_consensus_call",
    "ref",
    "cnv_excluded_regions.bed"
  ),
  col_names = c("chrom", "start", "end")
) %>%
  # Reformat the chromosome variable to drop the "chr"
  mutate(chrom = factor(gsub("chr", "", chrom),
                               levels = c(1:22, "X", "Y")
  )) %>%
  filter(
    # Drop CNVs that don't have chromosome labels
    !is.na(chrom),
    # Drop sex chromosomes
    !(chrom %in% c("X", "Y", "M"))
  )


# Set up uncallable regions as GenomicRanges. 
uncallable_ranges <- GenomicRanges::GRanges(
  seqnames = uncallable_bed$chrom,
  ranges = IRanges::IRanges(
    start = uncallable_bed$start,
    end = uncallable_bed$end
  )
)


## Call bin CN statuses for each sample

# Set up binned genome ranges. 

# Set up bins of ~1Mb size
bins <- GenomicRanges::tileGenome(
  chr_sizes_vector,
  tilewidth = 1e6
)
# Uncompress these ranges
bins <- unlist(bins)

# Get a vector of the biospecimen IDs
sample_ids <- unique(seg_data$ID)

# Read in the calculation file
bin_calls_df <- read_tsv(file.path(original_dir, "results", "cn_status_bp_per_bin.tsv"))


## Set up heatmap annotation objects

# Make color key. 
color_key <- structure(
  c(divergent_col_hex, "#9932CC", "#fed8b1"),
  names = c("loss", "neutral", "gain", "unstable", "uncallable")
)



### Make column annotation object
chrs <- paste0("chr", S4Vectors::decode(bins@seqnames))
chrs <- factor(chrs, levels = paste0("chr", 1:22))

# Make a key for assigning alternating colors to the chromosomes
chr_colors <- rep(c("white", "black"), length.out = length(unique(chrs)))
names(chr_colors) <- unique(chrs)

# Get coordinate start, end, midpoint positions
chr_start <- match(unique(chrs), chrs)
chr_end <- chr_start + summary(chrs)
mid_points <- floor((chr_start + chr_end) / 2)


# Make chromosomal labeling `HeatmapAnnotation` object.


# Make text labels for chromosome text
chr_text <- ComplexHeatmap::anno_mark(
  at = mid_points,
  labels = levels(chrs),
  which = "column",
  side = "bottom",
  labels_rot = 45,
  labels_gp = grid::gpar(cex = 0.65)
)

# Create the Heatmap annotation object
chr_annot <- ComplexHeatmap::HeatmapAnnotation(
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
histologies <- metadata %>%
  drop_na(cancer_group_display) %>%
  # No other ? TODO!
  filter(cancer_group_display != "Other") %>%
  # IDs of interest 
  filter(Kids_First_Biospecimen_ID %in% bin_calls_df$biospecimen_id) %>%
  # count, arrange, and factors for ordering
  count(cancer_group_display, name = "cancer_group_n") %>%
  arrange(-cancer_group_n) %>%
  # set an order column
  mutate(cancer_group_order = 1:n()) %>%
  inner_join(metadata) %>%
  # have to filter again after the join...
  filter(Kids_First_Biospecimen_ID %in% bin_calls_df$biospecimen_id) %>%
  # new column for display where N's are given
  mutate(cancer_group_display_n = glue::glue("{cancer_group_display} (N = {cancer_group_n})")) %>%
  # keep only columns we need moving forward
  select(Kids_First_Biospecimen_ID, 
         cancer_group_display_n, #labels
         cancer_group_order,  #order
         # color
         cancer_group_hex) %>%
 # Reorder display based on order
 mutate(cancer_group_display_n = forcats::fct_reorder(cancer_group_display_n, cancer_group_order)) %>%
 # Make sure they rows also in cancer_group_order
 arrange(cancer_group_display_n) %>%
 # Store as rownames
 column_to_rownames("Kids_First_Biospecimen_ID") %>%
 # We don't want this to actually be displayed on the heatmap though
 select(-cancer_group_order)

 
# Make a color key that's formatted for ComplexHeatmap
# Get a distinct version of the color keys
histologies_color_key_df <- histologies %>%
  dplyr::select(cancer_group_display_n, cancer_group_hex) %>%
  dplyr::distinct()

# Make color key specific to these samples
histologies_color_key <- histologies_color_key_df$cancer_group_hex
names(histologies_color_key) <- histologies_color_key_df$cancer_group_display_n

# Get coordinate start positions
hist_start <- match(names(histologies_color_key), histologies$cancer_group_display_n)

# Get coordinate end positions for each histology group
hist_end <- hist_start + summary(histologies$cancer_group_display_n)

# Get mid points of 
mid_points <- floor((hist_start + hist_end) /2)


# Make text labels for chromosome text
hist_text <- ComplexHeatmap::anno_mark(
  at = mid_points,
  labels = levels(histologies$cancer_group_display_n),
  which = "row",
  side = "right",
  labels_gp = grid::gpar(cex = 0.65),
  link_width = grid::unit(15, "mm")
)

# Create the Heatmap annotation object
hist_annot <- ComplexHeatmap::HeatmapAnnotation(
  df = data.frame(histologies),
  col = list(display_group = histologies_color_key),
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
bin_calls_mat <- bin_calls_mat[rownames(histologies), ]
if (all.equal(rownames(bin_calls_mat), rownames(histologies)) != TRUE) {
  stop("Bad data processing for Figure S2 Panel B")
}
 
## Assemble CN status heatmap
heatmap <- ComplexHeatmap::Heatmap(
  bin_calls_mat,
  name = "CN status",
  col = color_key,
  row_split = histologies$cancer_group_display_n,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  rect_gp = grid::gpar(col = "black", lwd = .0005),
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


# Print out heatmap. 
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")




























## Figure S3C-D ---------------------------
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


##### Figure S3C
fig_s3c <- merged_data %>%
  ggplot(aes(x = count_regions_any_conf_truncated, 
             y = cnv_breaks_count, 
             color = count_regions_any_conf_truncated)) +
  geom_jitter(width = 0.3, alpha = 0.9) +
  geom_boxplot(color = "black", alpha = 0, outlier.shape=NA) +
  ggsci::scale_color_simpsons() +
  theme(legend.position = "none") +
  xlab("Number of Chromothripsis Regions") + 
  ylab("Number of CNV Breaks") 
ggsave(figure_s3c_file, fig_s3c, width = 5, height = 3)

#### Figure S3D
fig_s3d <- merged_data %>%
  ggplot(aes(x = count_regions_any_conf_truncated, 
             y = sv_breaks_count, 
             color = count_regions_any_conf_truncated)) +
  geom_jitter(width = 0.3, alpha = 0.9) +
  geom_boxplot(color = "black", alpha = 0, outlier.shape=NA) +
  ggsci::scale_color_simpsons() +
  theme(legend.position = "none") +
  xlab("Number of Chromothripsis Regions") + 
  ylab("Number of SV Breaks") 
ggsave(figure_s3d_file, fig_s3d, width = 5, height = 3)


