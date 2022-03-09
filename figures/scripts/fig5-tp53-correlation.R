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

## Define output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig5", "correlation_plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Read in relevant data
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)
tp53_score <- readxl::read_excel(file.path(scores_dir, "TableS3-RNA-results-table.xlsx"), sheet = 1)
extend_score <- readxl::read_excel(file.path(scores_dir, "TableS3-RNA-results-table.xlsx"), sheet = 2)
breakpoint_den_data <- readr::read_tsv(file.path(breakpoint_den_dir, "cnv_breaks_densities.tsv")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = samples)

## Generate plots for tp53 scores vs. NormExtend Scores
tp53_extend <- left_join(tp53_score, extend_score)
tp53_extend$tp53_score <- as.numeric(tp53_extend$tp53_score)
tp53_extend$NormEXTENDScores_fpkm <- as.numeric(tp53_extend$NormEXTENDScores_fpkm)

### Save the plot
pdf(file.path(output_dir, "tp53_vs_extend_corr.pdf"))
p <- ggplot(data = tp53_extend , 
       aes(x = NormEXTENDScores_fpkm, y = tp53_score)) + 
  geom_point() + 
  xlab("NormEXTEND Scores") + 
  ylab("TP53 Scores") + 
  geom_smooth(method=lm, se=TRUE) + 
  stat_cor(method = "pearson", label.x = 0.05, label.y = 1, size = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
print(p)
dev.off()

## Generate plots for tp53 scores vs. breakpoint density
tp53_break <- left_join(tp53_score, breakpoint_den_data)
tp53_break$tp53_score <- as.numeric(tp53_break$tp53_score)
tp53_break$breaks_density <- as.numeric(tp53_break$breaks_density)

### Save the plot
pdf(file.path(output_dir, "tp53_vs_breakpoint_den_corr.pdf"))
p2 <- ggplot(data = tp53_break, 
       aes(x = breaks_density, y = tp53_score)) + 
  geom_point() + 
  xlab("Breakpoint Density") + 
  ylab("TP53 Scores") + 
  geom_smooth(method=lm, se=TRUE) + 
  stat_cor(method = "pearson", label.x = 0.05, label.y = 1, size = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16))
print(p2)
dev.off()
