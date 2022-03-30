# S. Spielman for CCDL 2022
#
# Makes pdf panels for supplementary Figure S2

library(tidyverse)


# Directories -------------------------------------------------------------------
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "supp", "figs2", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directores
analyses_dir <- file.path(root_dir, "analyses") 
snv_callers_dir <- file.path(analyses_dir, "snv-callers")
tmb_compare_dir <- file.path(analyses_dir, "tmb-compare")

# Scratch directory
scratch_dir <- file.path(root_dir, "scratch")

# Palette directory
palette_dir <- file.path(root_dir, "figures", "palettes")


## Define final PDF files ------------------------------------------------------
pbta_upset_pdf <- file.path(output_dir, "pbta_upset_plot.pdf")
tcga_upset_pdf <- file.path(output_dir, "tcga_upset_plot.pdf")
pbta_vaf_cor_matrix_pdf <- file.path(output_dir, "pbta_vaf_cor_matrix.pdf")
tcga_vaf_cor_matrix_pdf <- file.path(output_dir, "tcga_vaf_cor_matrix.pdf")
pbta_vaf_distribution_plot_pdf <- file.path(output_dir, "pbta_vaf_distribution_plot.pdf")
tcga_vaf_distribution_plot_pdf <- file.path(output_dir, "tcga_vaf_distribution_plot.pdf")
lancet_wxs_wgx_plot_pdf <- file.path(output_dir, "lancet_wxs_wgx_plot.pdf")

#A: pbta-vaf_cor_matrix.png
#B: pbta-vaf_distribution_plot.png
#C: pbta-upset-plot.png
#D: tcga-vaf-cor-matrix.png
#E: tcga-vaf-distribution-plot.png
#F: tcga-upset-plot.png
#G:http://htmlpreview.github.io/?https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/snv-callers/lancet-wxs-tests/lancet-paired-WXS-WGS.nb.html
#H/I are from tmb-compare


## Read in data and load items --------------------------------------------------------------

# Source function for the upset plots
source(file.path(snv_callers_dir, "util", "upset_plot.R"))

# Read in clinical data
histologies_df <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                           guess_max = 10000)

# Set up DB connections
con_pbta <- DBI::dbConnect(RSQLite::SQLite(), 
                           file.path(scratch_dir, "snv_db.sqlite"))
con_tcga <- DBI::dbConnect(RSQLite::SQLite(), 
                      file.path(scratch_dir, "tcga_snv_db.sqlite"))


# Define columns to join by
join_cols = c("Chromosome",
              "Start_Position",
              "Reference_Allele",
              "Allele",
              "Tumor_Sample_Barcode", 
              "Variant_Classification")


# Read in database tables, retaining only columns we need
strelka <- tbl(con_pbta, "strelka") %>% 
  select(join_cols, "VAF")

lancet <- tbl(con_pbta, "lancet") %>% 
  select(join_cols, "VAF")

mutect <- tbl(con_pbta, "mutect") %>% 
  select(join_cols, "VAF")

vardict <- tbl(con_pbta, "vardict") %>% 
  select(join_cols, "VAF")

# Source a script that will full join these and create `all_caller`
source(file.path(snv_callers_dir, "util", "full_join_callers.R"))


all_caller_df <- all_caller %>% 
  as.data.frame() %>% 
  rowid_to_column("index") %>%
  # Bring over VAF columns and the index
  select(starts_with("VAF_")) %>% 
  as.matrix()

# Store the indices as dimnames
dimnames(detect_mat)[[1]] <- all_caller_df$index

# Turn into detected or not
detect_mat <- !is.na(detect_mat)

# Plot from `detect_mat`
upset_png(detect_mat, 
          plot_file_path = pbta_upset_pdf,
          plot_file_type = "PDF")





