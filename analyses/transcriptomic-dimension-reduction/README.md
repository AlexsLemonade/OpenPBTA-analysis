# Dimension reduction: gene expression data

This module runs dimension reduction—namely PCA, UMAP, and t-SNE—on gene expression data from kallisto and RSEM.

To the run the whole pipeline and produce 4 plots

1. RSEM stranded data
2. RSEM poly-A data
3. kallisto stranded data
4. kallisto poly-A data

with panels for PCA, UMAP, and t-SNE that have points colored by broad histology, simply run:

using OPENPBTA_BASE_SUBTYPING=1 to run this module using the pbta-histologies-base.tsv from data folder while running molecular-subtyping modules for release.
```sh
OPENPBTA_BASE_SUBTYPING=1 dimension-reduction-plots.sh
```

OR by default uses pbta-histologies.tsv from data folder
```sh
bash dimension-reduction-plots.sh
```

## General usage

To document this module beyond what is in the individual files, we list some example applications or extensions below.

Running dimension reduction, creating lists of plots, and generating the multipanel plots can all be run separately using the appropriate Rscripts in `scripts/`:

```
generate-multipanel-plot.R
get-plot-list.R
run-dimension-reduction.R
```

### Repeat dimension reduction step with different seeds

##### Command line

If we wanted to run UMAP multiple times with different seed values on the kallisto poly-A data, but _skip t-SNE_, we could do the following:

```sh
seeds=(309 4642 2019 2994 5783)

for seed in "${seeds[@]}"
do
  Rscript --vanilla scripts/run-dimension-reduction.R \
    --expression ../../data/pbta-gene-expression-kallisto.polya.rds \
    --metadata ../../data/pbta-histologies.tsv \
    --filename_lead kallisto_polyA_${seed} \
    --output_directory results \
    --seed ${seed} \
    --skip_tsne
done

```

Unfortunately, this approach will also run PCA 5 times. 

##### R

It's also possible to repeat this UMAP step 5 times within R:

```R
# load the required custom functions
source("util/dimension-reduction-functions.R")

# read in kallisto poly-A data and transpose
expression_file <- "../../data/pbta-gene-expression-kallisto.polya.rds"

# Read in expression data
expression_data <- readr::read_rds(expression_file) %>%
  as.data.frame()

# the first column will be either transcript or column ids
# depending on the expression file
feature_identifier <- colnames(expression_data)[1]

# drop any columns that contain other identifers
expression_data <- expression_data %>%
  dplyr::select(!!rlang::sym(feature_identifier), dplyr::starts_with("BS_")) %>%
  dplyr::filter(complete.cases(.)) %>%
  tibble::column_to_rownames(var = feature_identifier)

# Filter out low count genes
genes_to_keep <- rowSums(expression_data) >= 100
expression_data <- expression_data[genes_to_keep, ]

# Transpose the data
transposed_exp_data <- t(expression_data)

seeds <- c(309, 4642, 2019, 2994, 5783)
for (seed in seeds) {
  perform_dimension_reduction(transposed_expression_matrix = transposed_exp_data,
                              method = "UMAP",
                              seed = seed,
                              model_filename = paste("kallisto_polyA", seed, "UMAP.RDS"),
                              output_directory = "results",
                              neighbors_parameter = 15)
}
```

### Generate plots colored by short histology

To add multipanel plots with points colored by `short_histology` (must be a column in `pbta-histologies.tsv`), we can run the following:

```sh
# Change the point color variable when
COLOR=short_histology bash 02-get-dimension-reduction-plot-lists.sh

# This will make plots for anything with data
bash 03-multipanel-plots.sh
```

This will plot the first two variables for each approach.
Currently, this does not matter for the way t-SNE and UMAP are implemented.
The PCA plots are limited to PC1 and PC2.

### Plot additional principal components

If we wanted to plot additional PCs for the RSEM stranded data in R, we can take advantage of the custom `plot_dimension_reduction` function and use the following R code:

```R
# load the required custom functions
source("util/dimension-reduction-functions.R")

# read in the data.frame that has the PCA info and the metadata
pca_df <- readr::read_tsv("results/rsem_stranded_pca_scores_aligned.tsv")

# plot itself - PC7 and PC10
plot_dimension_reduction(aligned_scores_df = pca_df,
                         point_color = "broad_histology",
                         x_label = "PC7",
                         y_label = "PC10",
                         score1 = 7,
                         score2 = 10)
```

The `point_color` argument of `plot_dimension_reduction` can be changed to any variable that is a column in the original TSV file.

### Explore potential batch effects from the sequencing center 

The notebook `04-explore-sequencing-center-effects.Rmd` performs exploratory analyses and visualization to roughly assess the extent to which sequencing center is expected to induce batch effects.