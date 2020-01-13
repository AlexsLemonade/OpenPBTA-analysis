## Immune Deconvolution

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to use the R package immunedeconv to quantify and compare various immune cell types in the tumor microenvironment (TME) across various PBTA histologies.
The package `immunedeconv`, provides six deconvolution methods: xCell (n = 64; immune and non-immune cell types), CIBERSORT (n = 22; with and without absolute mode), TIMER (n = 6), EPIC (n = 6), quanTIseq (n = 10) and MCP-Counter (n = 8). 
We use two deconvolution methods: xCell with 64 cell types and CIBERSORT (abs.) with 22 cell types to get a more fine grained deconvolution of the tumor microenvironment.
CIBERSORT requires two files that are available upon request from https://cibersort.stanford.edu/.
We recommend using MCP-Counter instead of CIBERSORT (abs.) when these files are not available to the user.

### Choice of method:

We chose xCell as the method of choice because it: 
1) is able to deconvolute the maximum number of immune and non-immune cell types 
2) is highly robust against background predictions and 
3) can reliably identify the presence of immune cells at low abundances (0-1% infiltration depending on the immune cell type).
To compare the immune scores predicted by xCell with that of other methods, we chose either CIBERSORT (abs.) and MCP-Counter as those two are the only methods in addition to xCell that output immune scores as arbitrary scores representing cell type abundance.
All three methods allow comparison between samples (inter-sample comparisons), between cell types (intra-sample comparisons) as well as between different cancer types (inter-histology comparisons).
We chose CIBERSORT (abs.) over MCP-Counter as it has the greatest number of overlapping immune cell types with that of xCell:

| method | overlapping_cell_type | overlap |
|--------|---------------|---------|
| cibersort_abs, mcp_counter, xcell | Monocyte, Neutrophil, T cell CD8+ | 3 |
| mcp_counter, xcell | B cell, Cancer associated fibroblast, Endothelial cell, Monocyte, Myeloid dendritic cell, Neutrophil, NK cell, T cell CD8+ | 8 |
| cibersort_abs, xcell | B cell memory, B cell naive, B cell plasma, Eosinophil, Macrophage M1, Macrophage M2, Monocyte, Myeloid dendritic cell activated, Neutrophil, T cell CD4+ naive, T cell CD8+, T cell gamma delta, T cell regulatory (Tregs) | 13 |

### Analysis scripts

#### 01-immune-deconv.R

1. Inputs from data download

```
pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
```

2. Function

This script deconvolutes immune cell types using xCell and another method, which in our case is CIBERSORT (abs.).
The other recommended method is MCP-Counter (see above for explanation). 

3. Output: 

`results/deconv-output.RData`

The results in the RData object are predicted immune scores per cell type per input sample.
The predicted scores are not actual cell fractions but arbitrary scores which can be compared within samples, across samples and/or between various cancer types.
Depending on the requirement, a user could use these to create various plots like scatterplot, corrplot, barplot, boxplot etc. 

#### 02-summary-plots.R 

1. Input

`results/deconv-output.RData`

2. Function:

This script creates summary correlation plots and heatmaps from predicted immune scores.

3. Output

* `plots/corrplot_xCell_vs_CIBERSORT(abs.).png`: a correlation plot of predicted immune scores across 13 common cell types from xCell and CIBERSORT (abs.) 
* `plots/heatmap_xCell.png`: heatmap of average immune scores per cell type per histology for xCell.
* `plots/heatmap_CIBERSORT(abs.).png`: heatmap of average immune scores per cell type per histology for the second method, CIBERSORT in our case.

### Running the analysis

```sh
bash run-immune-deconv.sh
```



