# Run Jin for D3b 
# Generate correlation plots of TP53 vs. NormEXTEND and breakpoint density

library(tidyverse)
library(readxl)
library(ggpubr)

## Define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analyses_dir <- file.path(root_dir, "analyses")
scores_dir <- file.path(root_dir, "tables", "results")
breakpoint_den_dir <- file.path(analyses_dir, "chromosomal-instability", "breakpoint-data")
cg_display_dir <- file.path(root_dir, "figures", "palettes")

## Define output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig5", "correlation_plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Read in relevant data
# read in histologies file
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)
# read in cancer group display file
palette_info <- readr::read_tsv(file.path(cg_display_dir, "broad_histology_cancer_group_palette.tsv")) %>%
  dplyr::select(broad_histology, cancer_group, cancer_group_display)
# join the information about palette
histologies_df_palette <- left_join(histologies_df,  palette_info)

tp53_score <- readxl::read_excel(file.path(scores_dir, "TableS3-RNA-results-table.xlsx"), sheet = 1)
extend_score <- readxl::read_excel(file.path(scores_dir, "TableS3-RNA-results-table.xlsx"), sheet = 2)
breakpoint_den_data <- readr::read_tsv(file.path(breakpoint_den_dir, "cnv_breaks_densities.tsv")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = samples)

## Generate plots for tp53 scores vs. NormExtend Scores
tp53_extend <- left_join(tp53_score, extend_score)
tp53_extend$tp53_score <- as.numeric(tp53_extend$tp53_score)
tp53_extend$NormEXTENDScores_fpkm <- as.numeric(tp53_extend$NormEXTENDScores_fpkm)

# add cancer group to the dataframe
tp53_extend <- tp53_extend %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_RNA, tp53_score, NormEXTENDScores_fpkm) %>%
  dplyr::left_join(histologies_df_palette %>% 
                     dplyr::select(sample_id, cancer_group_display)) %>%
  distinct()

### Save the plot
pdf(file.path(output_dir, "tp53_vs_extend_corr.pdf"))
p <- ggplot(data = tp53_extend , 
       aes(x = NormEXTENDScores_fpkm, y = tp53_score),
       ggtheme = theme_pubr()) + 
  geom_point() + 
  xlab("NormEXTEND Scores") + 
  ylab("TP53 Scores") + 
  geom_smooth(method=lm, se=TRUE) + 
  stat_cor(method = "pearson", label.x = 0.05, label.y = 1, size = 2.5) + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 5)) +
  facet_wrap(~cancer_group_display)
print(p)
dev.off()

## Generate plots for tp53 scores vs. breakpoint density
tp53_break <- left_join(tp53_score, breakpoint_den_data)
tp53_break$tp53_score <- as.numeric(tp53_break$tp53_score)
tp53_break$breaks_density <- as.numeric(tp53_break$breaks_density)

# add cancer group to the dataframe
tp53_break <- tp53_break %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_DNA, tp53_score, breaks_density) %>%
  dplyr::left_join(histologies_df_palette %>% 
                     dplyr::select(sample_id, cancer_group_display)) %>%
  distinct()

### Save the plot
pdf(file.path(output_dir, "tp53_vs_breakpoint_den_corr.pdf"))
p2 <- ggplot(data = tp53_break, 
       aes(x = breaks_density, y = tp53_score)) + 
  geom_point() + 
  xlab("Breakpoint Density") + 
  ylab("TP53 Scores") + 
  geom_smooth(method=lm, se=TRUE) + 
  stat_cor(method = "pearson", label.x = 0.01, label.y = 1, size = 2.5) + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 5)) +
  facet_wrap(~cancer_group_display)
print(p2)
dev.off()
