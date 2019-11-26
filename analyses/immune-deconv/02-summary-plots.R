# Author: Komal S. Rathi
# Date: 11/11/2019
# Function: Summarise results and create plots

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(corrplot))

# source plotting theme
source('analyses/immune-deconv/pubTheme.R')

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Immunedeconv output from 01-immune.deconv.R (.RData)"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output PDF (.pdf)")
)

# Example Run:
# Rscript analyses/immune-deconv/02-summary-plots.R \
# -i 'analyses/immune-deconv/deconv-output.RData' \
# -o 'analyses/immune-deconv/deconv-summary.pdf'

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
deconvout <- opt$input
output <- opt$output
load(deconvout) 

# first, define a function to create heatmap of
# average immune scores per histology per cell type
create.heatmap <- function(deconv.method, title) {
  
  # create labels: count of samples per histology
  annot <- deconv.method %>%
    select(broad_histology, sample) %>%
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

# next, plot a correlation heatmap between xCell and Cibersort
# only take common cell types between both methods
common.types <- intersect(cibersort_abs$cell_type, xcell$cell_type) 
cibersort_abs.sub <- cibersort_abs %>%
  filter(cell_type %in% common.types) %>%
  mutate(cibersort = fraction) %>%
  select(-c(method, fraction))
xcell.sub <- xcell %>%
  filter(cell_type %in% common.types) %>%
  mutate(xcell = fraction) %>%
  select(-c(method, fraction))
total <- merge(xcell.sub, cibersort_abs.sub, by = c("sample","cell_type", "broad_histology"))

# Overall correlation: 0.12
print(paste("Overall correlation: ", round(cor(total$xcell, total$cibersort), 2))) 

# labels
total.labels <- total %>%
  select(broad_histology, sample) %>%
  unique() %>%
  group_by(broad_histology) %>%
  dplyr::summarise(label = n()) %>%
  mutate(label = paste0(broad_histology,' (',label,')'))

# add labels to actual data
total <- merge(total, total.labels, by = 'broad_histology')

# calculate correlation per cell type per histology
total <- total %>% 
  group_by(cell_type, label) %>%
  dplyr::summarise(corr = cor(xcell, cibersort)) %>%
  spread(key = label, value = corr) %>% 
  column_to_rownames('cell_type') %>%
  replace(is.na(.), 0)

# create correlation heatmap
pdf(file = output, onefile = TRUE, width = 13, height = 8)
corrplot(t(total), method = "circle", type = 'full', win.asp = 0.5, 
         addCoef.col = "black", number.cex = .5,
         is.corr = FALSE, tl.cex = 0.8, mar = c(0, 0, 0, 5), 
         title = "\n\n\n\nCorrelation matrix (xCell vs Cibersort)")

# lastly, create heatmaps for both deconvolution methods
# add to the same file as above
create.heatmap(deconv.method = xcell, title = 'xCell')
create.heatmap(deconv.method = cibersort_abs, title = 'Cibersort')
dev.off()

