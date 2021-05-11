# Author: Komal S. Rathi
# Function: Summarise results and create plots

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))

# source plotting theme and heatmap functions
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "immune-deconv", "util", "pubTheme.R"))
source(file.path(root_dir, "analyses", "immune-deconv", "util", "heatmap_by_histology.R"))
source(file.path(root_dir, "analyses", "immune-deconv", "util", "heatmap_by_molecular_subtype.R"))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Immunedeconv output from 01-immune.deconv.R (.RData)"),
  make_option(c("-o",  "--output_dir"), type = "character", help = "Output directory")
)

# Example Run:
# Rscript 02-summary-plots.R \
# -i 'results/deconv-output.RData' \
# -o 'plots

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
deconvout <- opt$input
output_dir <- opt$output_dir
load(deconvout) 

# deconvolution method
deconv_method <- unique(deconv_output$method)

# create heatmap of average immune scores per cell type per histology
output_file <- file.path(output_dir, paste0("heatmap_", deconv_method, "_by_histology.pdf"))
heatmap_by_histology(deconv_output = deconv_output, output_file = output_file)

# create heatmap of average immune scores per cell type per molecular subtype per histology 
output_file <- file.path(output_dir, paste0("heatmap_", deconv_method, "_by_molecular_subtype.pdf"))
pdf(output_file, width = 15, height = 8)
plyr::d_ply(deconv_output, 
            .variables = "broad_histology", 
            .fun = function(x) heatmap_by_molecular_subtype(deconv_output = x))
dev.off()
