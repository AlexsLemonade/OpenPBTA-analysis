# Run Jin and Jo Lynne Rokita 
#
# Generate figures with UMAP results 

### Load libraries
library(ggplot2)
library(tidyverse)

### Define directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", "release-v21-20210820")
figure_script_dir <- file.path(root_dir, "figures", "scripts")
dim_red_dir <- file.path(root_dir, "analyses", "transcriptomic-dimension-reduction")

# define output directory
plots_dir <- file.path(root_dir, "figures", "pdfs", "supp")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# source the function for generating plots
source(file.path(dim_red_dir, "util", "dimension-reduction-functions.R"))

# UMAP data file
dat <- readRDS(file.path(dim_red_dir, 
                         "plots", "plot_data",
                         "rsem_stranded_log_broad_histology_multiplot_list.RDS"))

# histology file
histology_df <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"), 
                                guess_max = 10000)

# extract UMAP results
umap_results <- dat$UMAP$data %>%
  dplyr::select(Kids_First_Biospecimen_ID, X1, X2)
umap_df <- umap_results %>%
  dplyr::inner_join(histology_df) %>%
  select(Kids_First_Biospecimen_ID, X1, X2, broad_histology, cancer_group, molecular_subtype)

# possible colors for the palette
palette_OkabeIto <- c("#E69F00", 
                      "#56B4E9", 
                      "#009E73", 
                      "#F0E442", 
                      "#0072B2", 
                      "#D55E00", 
                      "#CC79A7", 
                      "#999999")

### Plot for HGG
# for HGG and DMG that has molecular subypte, keep the molecular subtype 
# for anything else, recode them as `Other CNS Tumor`
umap_df_hgg <- umap_df %>%
  dplyr::mutate(subtypes_to_plot = case_when(
    grepl("HGG", molecular_subtype) ~ gsub("HGG, ", "", molecular_subtype),
    grepl("DMG", molecular_subtype) ~ gsub("DMG, ", "", molecular_subtype),
    TRUE ~ "Other CNS tumor"
  )) %>%
  separate(subtypes_to_plot, c("hgat_subtypes", "tp53_status"), sep = ", ") %>% 
  dplyr::mutate(tp53_status = case_when(
    is.na(tp53_status) ~ "TP53 unchanged",
    TRUE ~ tp53_status
  ))

# reorder the levels for plotting
umap_df_hgg$tp53_status <- factor(umap_df_hgg$tp53_status, levels = c("TP53 unchanged",
                                                                      "TP53 activated",
                                                                      "TP53 loss"))

umap_df_hgg$hgat_subtypes <- factor(umap_df_hgg$hgat_subtypes, levels = c("H3 wildtype",
                                                                          "H3 K28",
                                                                          "H3 G35",
                                                                          "IDH",
                                                                          "To be classified",
                                                                          "Other CNS tumor"))

# save the figure
p <- plot_dimension_reduction(umap_df_hgg,
                              point_color = "hgat_subtypes",
                              point_shape = "tp53_status",
                              x_label = "Dimension 1",
                              y_label = "Dimension 2",
                              alpha_value = 0.7,
                              score1 = 2,
                              score2 = 3,
                              color_palette = c(sample(palette_OkabeIto, 4), "#656565", "#a9a9a9"))
# save the figures
pdf(file.path(plots_dir, "supp_umap_hgg.pdf"))
print(p)
dev.off()

### Plot for LGAT
# for LGG, we consider LGG, GNG and GNT all as LGG and we remove the prefixes to lumps groups together
# for anything else, recode them as `Other CNS Tumor`
umap_df_lgg <- umap_df %>%
  dplyr::mutate(subtypes_to_plot = case_when(
    grepl("LGG", molecular_subtype) ~ gsub("LGG, ", "", molecular_subtype),
    grepl("GNT", molecular_subtype) ~ gsub("GNT, ", "", molecular_subtype),
    grepl("GNG", molecular_subtype) ~ gsub("GNG, ", "", molecular_subtype),
    TRUE ~ "Other CNS tumor"
  )) %>% 
  dplyr::mutate(cdkn_status = case_when(
    grepl("CDKN2A/B", molecular_subtype) ~ "altered",
    TRUE ~ "not altered"
  )) %>% 
  dplyr::mutate(subtypes_to_plot = case_when(
    grepl("CDKN2A/B", subtypes_to_plot) ~ gsub(", CDKN2A/B", "", subtypes_to_plot),
    TRUE ~ subtypes_to_plot
  )) %>% 
  dplyr::mutate(lgat_subtypes = case_when(
    grepl("-somatic", subtypes_to_plot) & !grepl("-germline", subtypes_to_plot) ~ gsub("-somatic", "", subtypes_to_plot),
    grepl("-germline", subtypes_to_plot) & !grepl("-somatic", subtypes_to_plot) ~ gsub("-germline", "", subtypes_to_plot),
    grepl("-somatic, NF1-germline", subtypes_to_plot) ~ gsub("-somatic, NF1-germline", "", subtypes_to_plot),
    TRUE ~ subtypes_to_plot
  )) %>% 
  group_by(lgat_subtypes) %>% 
  dplyr::mutate(sample_n = n()) %>% 
  dplyr::mutate(colored_subtype = case_when(
    sample_n < 10 ~ "Other LGAT subtypes", 
    TRUE ~ lgat_subtypes
  )) %>% 
  dplyr::mutate(labels_subtype = case_when(
    sample_n < 10 ~ lgat_subtypes,
    TRUE ~ NA_character_
  ))

# reorder the levels for plotting
umap_df_lgg$cdkn_status <- factor(umap_df_lgg$cdkn_status, levels = c("not altered",
                                                                      "altered"))

umap_df_lgg$colored_subtype <- factor(umap_df_lgg$colored_subtype, levels = c("BRAF V600E",
                                                                              "KIAA1549-BRAF",
                                                                              "NF1",
                                                                              "other MAPK",
                                                                              "RTK",
                                                                              "wildtype",
                                                                              "Other LGAT subtypes",
                                                                              "To be classified",
                                                                              "Other CNS tumor"))
# generate a plot
p <- ggplot(umap_df_lgg, aes(x = X1,
                             y = X2,
                             color = colored_subtype,
                             shape = cdkn_status)) +
  ggplot2::scale_color_manual(values = c(sample(palette_OkabeIto, 6), "#000000", "#656565", "#a9a9a9")) + 
  ggplot2::geom_point(alpha = 0.7) +
  ggplot2::labs(x = "Dimension 1", y = "Dimension 2") +
  ggplot2::theme_bw() + 
  ggrepel::geom_text_repel(aes(label = labels_subtype),
                           size = 2)

# save the figure
pdf(file.path(plots_dir, "supp_umap_lgg.pdf"))
print(p)
dev.off()

### Plot for MB 
# for MB, we keep subtypes and recode everything else as `Other CNS Tumor`
umap_df_mb <- umap_df %>%
  dplyr::mutate(mb_subtypes = case_when(
    grepl("MB", molecular_subtype) ~ gsub("MB, ", "", molecular_subtype),
    TRUE ~ "Other CNS tumor"
  ))

# reorder the levels for plotting
umap_df_mb$mb_subtypes <- factor(umap_df_mb$mb_subtypes, levels = c("Group3",
                                                                    "Group4",
                                                                    "SHH",
                                                                    "WNT",
                                                                    "To be classified",
                                                                    "Other CNS tumor"))

# save the plot
p <- plot_dimension_reduction(umap_df_mb,
                              point_color = "mb_subtypes",
                              x_label = "Dimension 1",
                              y_label = "Dimension 2",
                              alpha_value = 0.7,
                              score1 = 2,
                              score2 = 3,
                              color_palette = c(sample(palette_OkabeIto, 4), "#656565", "#a9a9a9"))
# save the figure
pdf(file.path(plots_dir, "supp_umap_mb.pdf"))
print(p)
dev.off()

### Plot for EPN
# for EPN, we keep subtypes and recode everything else as `Other CNS Tumor`
umap_df_epn <- umap_df %>%
  dplyr::mutate(epn_subtypes = case_when(
    grepl("EPN", molecular_subtype) ~ gsub("EPN, ", "", molecular_subtype),
    TRUE ~ "Other CNS tumor"
  ))

# reorder the levels for plotting
umap_df_epn$epn_subtypes <- factor(umap_df_epn$epn_subtypes, levels = c("ST RELA",
                                                                        "ST YAP1",
                                                                        "PF A",
                                                                        "H3 K28",
                                                                        "To be classified",
                                                                        "Other CNS tumor"))
# save the plots and output in the console
p <- plot_dimension_reduction(umap_df_epn,
                              point_color = "epn_subtypes",
                              x_label = "Dimension 1",
                              y_label = "Dimension 2",
                              alpha_value = 0.7,
                              score1 = 2,
                              score2 = 3,
                              color_palette = c(sample(palette_OkabeIto, 4), "#656565", "#a9a9a9"))
# save the figure
pdf(file.path(plots_dir, "supp_umap_epn.pdf"))
print(p)
dev.off()


