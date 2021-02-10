# Author: Komal S. Rathi
# Date: 11/11/2019
# Function: Summarise results and create plots

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(corrplot))

# source plotting theme
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "immune-deconv",
                 "util", "pubTheme.R"))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Immunedeconv output from 01-immune.deconv.R (.RData)"),
  make_option(c("-o",  "--output"), type = "character", help = "Output directory")
)

# Example Run:
# Rscript analyses/immune-deconv/02-summary-plots.R \
# -i 'analyses/immune-deconv/results/deconv-output.RData' \
# -o 'analyses/immune-deconv/plots

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
deconvout <- opt$input
output <- opt$output
load(deconvout) 

# add molecular subtype info
deconv.res$broad_histology <- ifelse(is.na(deconv.res$molecular_subtype),
                                     deconv.res$broad_histology,
                                     paste0(deconv.res$broad_histology, '-', deconv.res$molecular_subtype))


# extract names of the methods used
methods <- unique(deconv.res$method)
method1.name <- methods[1]
method2.name <- methods[2]

# split results of method 1 and 2
method1 <- deconv.res %>%
  filter(method == method1.name)
method2 <- deconv.res %>%
  filter(method == method2.name)

# first, define a function to create heatmap of
# average immune scores per histology per cell type
create.heatmap <- function(deconv.method, title, fileout) {

  # assign labels
  non.brain.tumors <- c("Histiocytic tumor", "Lymphomas")

  # create labels: count of samples per histology
  annot <- deconv.method %>%
    dplyr::select(broad_histology, sample) %>%
    unique() %>%
    group_by(broad_histology) %>%
    summarise(label = n()) %>%
    mutate(label = paste0(broad_histology,' (',label,')'))

  # add labels to actual data
  deconv.method <- merge(deconv.method, annot, by = 'broad_histology')

  # calculate average scores per cell type per histology
  deconv.method <- deconv.method %>%
    filter(!cell_type %in% c("microenvironment score", "stroma score", "immune score")) %>%
    group_by(cell_type, label) %>%
    dplyr::summarise(mean = mean(fraction)) %>%
    # convert into matrix of cell type vs histology
    spread(key = label, value = mean) %>%
    column_to_rownames('cell_type')

  # plot non-brain and brain tumors separately
  pdf(file = fileout, width = 15, height = 15)
  # non-brain tumors
  mat <- deconv.method %>%
    dplyr::select(grep(paste0(non.brain.tumors, collapse="|"), colnames(deconv.method), value = TRUE))
  if(ncol(mat) > 1){
    mat <- mat %>%
      rownames_to_column('celltype') %>%
      filter_at(vars(-celltype), any_vars(. != 0)) %>%
      column_to_rownames('celltype') %>%
      t() %>%
      pheatmap(fontsize = 10,
               scale = "column", angle_col = 45,
               main = "Average immune scores normalized by rows\nNon-Brain Tumors",
               annotation_legend = T, cellwidth = 15, cellheight = 15)
  }

  # brain tumors
  mat <- deconv.method %>%
    dplyr::select(grep(paste0(non.brain.tumors, collapse="|"), colnames(deconv.method), invert = TRUE, value = TRUE))
  if(ncol(mat) > 1){
    mat <- mat %>%
      rownames_to_column('celltype') %>%
      filter_at(vars(-celltype), any_vars(. != 0)) %>%
      column_to_rownames('celltype') %>%
      t() %>%
      pheatmap(fontsize = 10,
               scale = "column", angle_col = 45,
               main = "Average immune scores normalized by rows\nBrain Tumors",
               annotation_legend = T, cellwidth = 15, cellheight = 15)
  }

  dev.off()
}

# next, plot a correlation heatmap between xCell and the second specified method
# only take common cell types between both methods
common.types <- intersect(method1$cell_type, method2$cell_type)
method1.sub <- method1 %>%
  filter(cell_type %in% common.types) %>%
  mutate(!!method1.name := fraction) %>%
  dplyr::select(-c(method, fraction))
method2.sub <- method2 %>%
  filter(cell_type %in% common.types) %>%
  mutate(!!method2.name := fraction) %>%
  dplyr::select(-c(method, fraction))
total <- merge(method1.sub, method2.sub, by = c("sample","cell_type", "broad_histology", "molecular_subtype"))

# Overall correlation
avg.cor <- round(cor(total[,method1.name], total[,method2.name]), 2)
print(paste("Overall Pearson Correlation: ", avg.cor))

# labels
total.labels <- total %>%
  dplyr::select(broad_histology, sample) %>%
  unique() %>%
  group_by(broad_histology) %>%
  dplyr::summarise(label = n()) %>%
  mutate(label = paste0(broad_histology,' (',label,')'))

# add labels to actual data
total <- merge(total, total.labels, by = 'broad_histology')

# calculate correlation per cell type per histology
total <- total %>%
  group_by(cell_type, label) %>%
  dplyr::summarise(corr = cor(!!sym(method1.name), !!sym(method2.name))) %>%
  spread(key = label, value = corr) %>%
  column_to_rownames('cell_type') %>%
  replace(is.na(.), 0)

# replace space from method names for output filename
m1 <- gsub(" ","",method1.name)
m2 <- gsub(" ","",method2.name)

# create correlation plot for overlapping cell types between both methods
pdf(file = file.path(output, paste0("corrplot_", m1, "_vs_", m2, ".pdf")),
    width = 16, height = 20)
corrplot(t(total), method = "circle", type = 'full', win.asp = 0.5,
         addCoef.col = "#888888", number.cex = .7,
         tl.col = "black", number.font = 2,
         is.corr = FALSE, tl.cex = 0.8,
         mar = c(0, 0, 0, 5),
         cl.cex = .8,
         title = paste0("\n\n\n\nCorrelation matrix (",
                        method1.name, " vs ", method2.name, ")\n",
                        "Overall Pearson Correlation: ", avg.cor))
dev.off()

# create heatmaps of average immune scores per histology per cell type
# method1
create.heatmap(deconv.method = method1, title = method1.name,
               fileout = file.path(output, paste0("heatmap_", m1, ".pdf")))

# method2
create.heatmap(deconv.method = method2, title = method2.name,
               fileout = file.path(output, paste0("heatmap_", m2, ".pdf")))
