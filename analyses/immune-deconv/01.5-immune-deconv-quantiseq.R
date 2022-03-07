# This script will run QuantiSEQ via the immunedeconv package.
# We will perform with BOTH stranded and polya separately and then compare. If same, we pool.

# Find the root directory of this repository
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(immunedeconv))

option_list <- list(
  make_option(c("-p", "--polyaexprs"), type = "character",
              help = "PolyA Expression data: HUGO symbol x Sample identifiers (.RDS)"),
  make_option(c("-s", "--strandedexprs"), type = "character",
              help = "Stranded Expression data: HUGO symbol x Sample identifiers (.RDS)"),
  make_option(c("-c", "--clin"), type = "character",
              help = "Clinical file (.TSV)"),
  make_option(c("-m", "--method"), type = "character",
              default = "xcell",
              help = "Deconvolution Method"),
  make_option(c("-o","--outputfile"), type = "character",
              help = "Deconv Output (.RData)")
)

# Example Run:
# Rscript 01-immune-deconv.R \
# --polyaexprs '../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds' \
# --strandedexprs '../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds' \
# --clin '../../data/pbta-histologies.tsv' \
# --method 'xcell' \
# --outputfile 'results/deconv-output.RData'

# parse parameters
#opt <- parse_args(OptionParser(option_list = option_list))
polya <- '../../data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds' #opt$polyaexprs
stranded <- '../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds' #opt$strandedexprs
clin_file <- '../../data/pbta-histologies.tsv' #opt$clin
deconv_method <- "quantiseq" #opt$method
#output_file <- #opt$outputfile

#### Check model parameter - must be in deconvolution_methods (immunedeconv accepted options)
if (!is.null(deconv_method)){
  if (!(deconv_method %in% deconvolution_methods)) {
    stop( paste(c("Specified method not available. Must be one of the following: ", deconvolution_methods), collapse=" ") )
  }
}

# merge expression from polya and stranded data on common genes
polya <- readRDS(polya)
stranded <- readRDS(stranded)

# read clinical data
clin_file <- readr::read_tsv(clin_file, guess_max = 10000)

# Import standard color palettes for project
histology_label_mapping <- readr::read_tsv(
file.path(root_dir, "figures", "palettes", "histology_label_color_table.tsv")) %>%
# Select just the columns we will need for plotting
dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes) %>%
# Reorder display_group based on display_order
dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))

# subset clinical
clin_file_sub  <- clin_file %>%
#filter(Kids_First_Biospecimen_ID %in% colnames(expr_input)) %>%
dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID") %>%
dplyr::select(Kids_First_Biospecimen_ID, broad_histology, display_group, molecular_subtype)

# deconvolute polya and stranded separately. 
# with 4 GB RAM, polya runs in 3-6 minutes, and stranded in 20 minutes or so.
# `deconvolute` function args we need:
## Default is `tumor = TRUE`, which is what we want
## Default is `scale_mRNA = TRUE`. This performs correction of cell-type-specific mRNA content bias.
## we want to keep scale_mRNA = TRUE since we are NOT working with simulated data
result_polya <- deconvolute(gene_expression = as.matrix(polya), method = "quantiseq")
result_stranded <- deconvolute(gene_expression = as.matrix(stranded), method = "quantiseq")

result_polya %>%
  gather(-cell_type, key = "sample", value = "score") %>%
  mutate(type = "polya") -> polya_wide
result_stranded %>%
  gather(-cell_type, key = "sample", value = "score") %>%
  mutate(type = "stranded") %>%
  bind_rows(polya_wide) -> results_quantiseq

write_tsv(results_quantiseq, "results/quantiseq_results.tsv")



load('results/deconv-output.RData')
results_xcell <- deconv_output %>% as_tibble()
results_xcell %>%
  filter(cell_type %in% unique(results_quantiseq$cell_type)) %>%
  select(-broad_histology, -display_group, -molecular_subtype, -method) %>%
  left_join(results_quantiseq, by = c("cell_type", "sample")) %>%
  rename(xcell = fraction, 
         quantiseq = score) -> data_compare_xcell_quantiseq

data_compare_xcell_quantiseq %>%
  filter(xcell < 1) %>%
  ggplot() +
  aes(x = xcell, y = quantiseq, color = cell_type) + 
  geom_point(alpha = 0.25)

data_compare_xcell_quantiseq %>%
  group_by(cell_type, type) %>%
  dplyr::summarize(pearson_correlation = cor(xcell, quantiseq, method = "pearson"), 
                   rank_correlation = cor(xcell, quantiseq, method = "spearman")) %>%
  knitr::kable()

