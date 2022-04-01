# Script for creating the Telomerase activities figure

library(stringr)
library(gridBase)
library(gridGraphics)
library(optparse)
library(forcats) # for fct_reorder()
library(ggpubr) # For pvalue_stat_manual

################################ File paths and reading data for making figure ###################################################################

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

Telomerase_PolyaCounts <- file.path(
  root_dir, "analyses", "telomerase-activity-prediction",
  "results", "TelomeraseScores_PTBAPolya_counts.txt"
) ### Variable representing Telomerase activity-prediction for Polya_counts

Telomerase_StdCounts <- file.path(
  root_dir, "analyses", "telomerase-activity-prediction",
  "results", "TelomeraseScores_PTBAStranded_counts.txt"
) ### Variable representing Telomerase activity-prediction for Stranded_counts

Telomerase_PolyaFpkm <- file.path(
  root_dir, "analyses", "telomerase-activity-prediction",
  "results", "TelomeraseScores_PTBAPolya_FPKM.txt"
) ### Variable representing Telomerase activity-prediction for PolyA_FPKM data

Telomerase_StdFpkm <- file.path(
  root_dir, "analyses", "telomerase-activity-prediction",
  "results", "TelomeraseScores_PTBAStranded_FPKM.txt"
) ### Variable representing Telomerase activity-prediction for Stranded_FPKM data



Histologies <- file.path(root_dir, "data", "pbta-histologies.tsv") ### Variable representing clinical data



palette_dir <- file.path(root_dir, "figures", "palettes")

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "fig4", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

telomerase_pdf <- file.path(output_dir, "Telomerase_Activities.pdf")


supplementary_telomerase_png <- file.path(root_dir, "figures", "pngs", "SuppTelomerase_Activities.png")

# Get palette for cancer group
cancer_group_palette <- readr::read_tsv(
  file.path(palette_dir, "broad_histology_cancer_group_palette.tsv")
) %>%
  dplyr::select(cancer_group, cancer_group_hex) %>%
  # Remove NA values -- a cancer group hex value will be NA only if the
  # cancer group is NA
  dplyr::filter(complete.cases(.))

# Make color palette suitable for use with ggplot
annotation_colors <- cancer_group_palette$cancer_group_hex
names(annotation_colors) <- cancer_group_palette$cancer_group

# We need to map between biospecimen ID and cancer group
cancer_group_id_df <- readr::read_tsv(Histologies) %>%
  dplyr::filter(!is.na(cancer_group)) %>%
  dplyr::select(Kids_First_Biospecimen_ID, cancer_group) %>%
  ## Renaming "Kids_First_Biospecimen_ID" as SampleID for comparison purpose
  dplyr::rename("SampleID" = "Kids_First_Biospecimen_ID")

TMScores1 <- read.table(Telomerase_StdFpkm, sep = "	", head = T) ## Reading Stranded FPKM telomerase scores
colnames(TMScores1)[colnames(TMScores1) == "NormEXTENDScores"] <- "NormEXTENDScores_Stranded_FPKM"

PTBA_GE_Standard_Histology <- merge(cancer_group_id_df, TMScores1, by = "SampleID") ### Merging Clinical data with the Telomerase scores

TMScores2 <- read.table(Telomerase_PolyaFpkm, sep = "	", head = T)
colnames(TMScores2)[colnames(TMScores2) == "NormEXTENDScores"] <- "NormEXTENDScores_PolyA_FPKM"


TMScores3 <- read.table(Telomerase_StdCounts, sep = "	", head = T)
colnames(TMScores3)[colnames(TMScores3) == "NormEXTENDScores"] <- "NormEXTENDScores_StrandedCounts"

TMScores4 <- read.table(Telomerase_PolyaCounts, sep = "	", head = T)
colnames(TMScores4)[colnames(TMScores4) == "NormEXTENDScores"] <- "NormEXTENDScores_PolyACounts"

PBTA_PolyA_TMScores <- merge(TMScores2, TMScores4, by = "SampleID")
PBTA_Stranded_TMScores <- merge(TMScores1, TMScores3, by = "SampleID")

########################################## Figure A and B dataframe compilation #########################################################

# This data frame will be used for most histology plots
Stranded_Histology <- PTBA_GE_Standard_Histology

#########################################  Saving Figure in PNG format


## Figure for main text: Boxplots
theme_set(theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 60, size = 6, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 0),
    axis.title.y = element_text(size = 8),
    legend.position = "none",
    legend.key.size = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ))

P1 <- ggplot(Stranded_Histology , aes(
  x = fct_reorder(cancer_group, NormEXTENDScores_Stranded_FPKM, .desc = TRUE),
  y = NormEXTENDScores_Stranded_FPKM
)) +
  geom_boxplot(
    size = 0.2, notch = FALSE, outlier.size = 0, outlier.shape = NA,
    aes(color = cancer_group, fill = cancer_group), alpha = 0.65
  ) +
  geom_jitter(shape = 16, cex = 0.2, aes(color = cancer_group),
              alpha = 0.75) +
  scale_fill_manual(values = annotation_colors, aesthetics = c("colour", "fill")) +
  ylab("EXTEND Scores (Stranded FPKM)") +
  xlab("Cancer Group")

ggsave(plot = P1, telomerase_pdf, dpi = 1200, units = "in",
       width = 8, height = 4)

## Figure for SI: scatterplots
png(supplementary_telomerase_png, width = 4, height = 2, units = "in", res = 1200)

theme_set(theme_classic() +
  theme(
    plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(angle = 25, size = 6, vjust = 1, hjust = 1), axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 8),
    legend.position = "none",
    legend.key.size = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  ))


P1 <- ggscatter(PBTA_PolyA_TMScores,
  x = "NormEXTENDScores_PolyACounts", y = "NormEXTENDScores_PolyA_FPKM", color = "red", size = 0.2,
  add = "reg.line", # Add regression line
  add.params = list(color = "black", fill = "grey", size = 0.5), # Customize reg. line
  conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", size = 2)


P2 <- ggscatter(PBTA_Stranded_TMScores,
  x = "NormEXTENDScores_StrandedCounts", y = "NormEXTENDScores_Stranded_FPKM", color = "red", size = 0.2,
  add = "reg.line", # Add regression line
  add.params = list(color = "black", fill = "grey", size = 0.5), # Customize reg. line
  conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", size = 2)



grid.newpage()
# Create layout : nrow = 2, ncol =2
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
# A helper function to define a region on the layout
define_region <- function(row, col) {
  viewport(layout.pos.row = row, layout.pos.col = col)
}



print(ggpar(P1, font.xtickslab = 4, font.ytickslab = 4, font.x = 4, font.y = 4, xlab = "PolyA Counts EXTEND Scores", ylab = "PolyA FPKM EXTEND Scores", title = "A", font.title = 7), vp = define_region(row = 1, col = 1))
print(ggpar(P2, font.xtickslab = 4, font.ytickslab = 4, font.x = 4, font.y = 4, xlab = "Stranded Counts EXTEND Scores ", ylab = "Stranded FPKM EXTEND Scores", title = "B", font.title = 7), vp = define_region(row = 1, col = 2))

dev.off()
