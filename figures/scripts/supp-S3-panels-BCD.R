# Stephanie J. Spielman for CCDL, 2022

# This script creates three S3 panels (B, C, D)


## Load libraries ------------------------
library(tidyverse)
library(ComplexHeatmap)


## Define paths  -----------------------
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses")
supp_figures_dir <- file.path(root_dir, "figures", "pdfs", "supp")

## Define output files
figure_s3b_file <- file.path(supp_figures_dir, "supp_figS3-B.pdf")
figure_s3c_file <- file.path(supp_figures_dir, "supp_figS3-C.pdf")
figure_s3d_file <- file.path(supp_figures_dir, "supp_figS3-D.pdf")

## Figure S3B ----------------------------
# We copy this figure over from the analysis module because it is created from an entire notebook
# which is too large to duplicate here
original_figure_file <- file.path(analyses_dir, "cnv-chrom-plot", "plots", "cn_status_heatmap.pdf")
file.copy(original_figure_file, figure_s3b_file)


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


