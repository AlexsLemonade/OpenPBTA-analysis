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
create.heatmap <- function(deconv.method, title) {
  
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
    summarise(mean = mean(fraction)) %>%
    # convert into matrix of cell type vs histology
    spread(key = label, value = mean) %>% 
    column_to_rownames('cell_type')
  
  # remove rows with all zeros (not allowed because we are scaling by row)
  deconv.method <- deconv.method[apply(deconv.method, 1, function(x) !all(x==0)),]
  
  title <- paste0(title,"\nAverage immune scores normalized by rows")
  pheatmap(mat = t(deconv.method), fontsize = 10, 
           scale = "column", angle_col = 45,
           main = title, annotation_legend = T, cellwidth = 15, cellheight = 15)
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
total <- merge(method1.sub, method2.sub, by = c("sample","cell_type", "broad_histology"))

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
png(filename = file.path(output, paste0("corrplot_", m1, "_vs_", m2, ".png")), 
    width = 13, height = 8, units = "in", res = 300)
corrplot(t(total), method = "circle", type = 'full', win.asp = 0.5, 
         addCoef.col = "black", number.cex = .5,
         is.corr = FALSE, tl.cex = 0.8, mar = c(0, 0, 0, 5), 
         title = paste0("\n\n\n\nCorrelation matrix (", 
                        method1.name, " vs ", method2.name, ")\n",
                        "Overall Pearson Correlation: ", avg.cor))
dev.off()

# create heatmaps of average immune scores per histology per cell type
# method1
png(filename = file.path(output, paste0("heatmap_", m1, ".png")), width = 13, height = 8, units = "in", res = 300)
create.heatmap(deconv.method = method1, title = method1.name)
dev.off()

# method2
png(filename = file.path(output, paste0("heatmap_", m2, ".png")), width = 10, height = 8, units = "in", res = 300)
create.heatmap(deconv.method = method2, title = method2.name)
dev.off()

