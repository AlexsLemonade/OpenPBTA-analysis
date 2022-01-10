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
dim_red_dir <- file.path(root_dir, "analysis", "transcriptomic-dimension-reduction")

# define output directory
plots_dir <- file.path(root_dir, "figures", "pdfs", "supp")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# source the function for generating plots
source(file.path(dim_red_dir, "util", "dimension-reduction-functions.R"))

### read in files
# UMAP data file
dat <- readRDS(file.path(dim_red_dir, 
                         "plots", "plot_data",
                         "rsem_stranded_log_broad_histology_multiplot_list.RDS"))

# histology file
histology_df <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"), 
                                guess_max = 10000)

### get UMAP results 
# extract UMAP results
umap_results <- dat$UMAP$data %>%
  dplyr::select(Kids_First_Biospecimen_ID, X1, X2)
umap_df <- umap_results %>%
  dplyr::inner_join(histology_df) %>%
  select(Kids_First_Biospecimen_ID, X1, X2, broad_histology, cancer_group, molecular_subtype)


### Plot for HGG
# for HGG and DMG that has molecular subypte, keep the molecular subtype 
# for anything else, recode them as `Other CNS Tumor`
umap_df_hgg <- umap_df %>%
  dplyr::mutate(subtypes_to_plot = case_when(
    grepl("HGG", molecular_subtype) & !grepl("To be classified", molecular_subtype) ~ molecular_subtype,
    grepl("DMG", molecular_subtype) & !grepl("To be classified", molecular_subtype) ~ molecular_subtype,
    TRUE ~ "Other CNS tumor"
  ))

# save the figure
p <- plot_dimension_reduction(umap_df_hgg,
                              point_color = "mol_alt",
                              point_shape = "subtype_code",
                              x_label = "Dimension 1",
                              y_label = "Dimension 2",
                              alpha_value = 0.3,
                              color_palette = NULL)
# save the figures
pdf(file.path(plots_dir, "supp_umap_hgg.pdf"))
print(p)
dev.off()

### Plot for LGAT
# for LGG, we consider LGG, GNG and GNT all as LGG and we remove the prefixes to lumps groups together
# for anything else, recode them as `Other CNS Tumor`
umap_df_lgg <- umap_df %>%
  dplyr::mutate(subtypes_to_plot = case_when(
    grepl("LGG", molecular_subtype) & !grepl("To be classified", molecular_subtype) ~ gsub("LGG, ", "", molecular_subtype),
    grepl("GNT", molecular_subtype) & !grepl("To be classified", molecular_subtype) ~ gsub("GNT, ", "", molecular_subtype),
    grepl("GNG", molecular_subtype) & !grepl("To be classified", molecular_subtype) ~ gsub("GNG, ", "", molecular_subtype),
    TRUE ~ "Other CNS tumor"
  ))

# save the figure
p <- plot_dimension_reduction(umap_df_lgg,
                              point_color = subtypes_to_plot,
                              x_label = "Dimension 1",
                              y_label = "Dimension 2",
                              alpha_value = 0.3,
                              color_palette = NULL)
# save the figures
pdf(file.path(plots_dir, "supp_umap_lgg.pdf"))
print(p)
dev.off()


### Plot for MB 
# for MB, we keep subtypes that are not `To be classfied` and recode everything else as `Other CNS Tumor`
umap_df_mb <- umap_df %>%
  dplyr::mutate(subtypes_to_plot = case_when(
    grepl("MB", molecular_subtype) & !grepl("To be classified", molecular_subtype) ~ molecular_subtype,
    TRUE ~ "Other CNS tumor"
  ))

# save the plot
p <- plot_dimension_reduction(umap_df_mb,
                              point_color = subtypes_to_plot,
                              x_label = "Dimension 1",
                              y_label = "Dimension 2",
                              alpha_value = 0.3,
                              color_palette = NULL)
# save the figure
pdf(file.path(plots_dir, "supp_umap_mb.pdf"))
print(p)
dev.off()

### Plot for EPN
# for EPN, we keep subtypes that are not `To be classfied` and recode everything else as `Other CNS Tumor`
umap_df_epn <- umap_df %>%
  dplyr::mutate(subtypes_to_plot = case_when(
    grepl("EPN", molecular_subtype) & !grepl("To be classified", molecular_subtype) ~ molecular_subtype,
    TRUE ~ "Other CNS tumor"
  ))

# save the plots and output in the console
p <- plot_dimension_reduction(umap_df_epn,
                              point_color = subtypes_to_plot,
                              x_label = "Dimension 1",
                              y_label = "Dimension 2",
                              alpha_value = 0.3,
                              color_palette = NULL)
# save the figure
pdf(file.path(plots_dir, "supp_umap_epn.pdf"))
print(p)
dev.off()


