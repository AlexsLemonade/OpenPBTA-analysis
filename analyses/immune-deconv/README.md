## Immune Deconvolution

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to use the R package `immunedeconv` to quantify and compare various immune cell types in the tumor microenvironment (TME) across various PBTA histologies. 
The package `immunedeconv`, provides six deconvolution methods: xCell (n = 64; immune and non-immune cell types), CIBERSORT (n = 22; with and without absolute mode), TIMER (n = 6), EPIC (n = 6), quanTIseq (n = 10) and MCP-Counter (n = 8). 

### Method selection

We chose xCell as the method of choice because it: 
1) is the most comprehensive deconvolution method and is able to deconvolute the maximum number of immune and non-immune cell types 
2) is highly robust against background predictions and 
3) can reliably identify the presence of immune cells at low abundances (0-1% infiltration depending on the immune cell type).

xCell outputs immune scores as arbitrary scores that represent cell type abundance. 
This allows comparison between samples (inter-sample comparisons), between cell types (intra-sample comparisons) as well as between different cancer types (inter-histology comparisons).

### Analysis scripts

#### 01-immune-deconv.R

1. Inputs from data download

```
pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
```

2. Function

This script deconvolutes immune cell types using `xCell`.

3. Output: 

`results/deconv-output.RData`

The results in the RData object are predicted immune scores per cell type per input sample. 
These scores are not actual cell fractions but arbitrary scores which can be compared within samples, across samples and/or between various cancer types. 
Depending on the user requirement, the output can be used to create various visualizations. 

#### 02-summary-plots.R 

1. Input

`results/deconv-output.RData`

2. Function:

This script creates heatmaps from predicted immune scores.

3. Output

* `plots/heatmap_xCell_by_histology.pdf`: heatmap of average immune scores per cell type per histology stratified by brain and non-brain tumors.
* `plots/heatmap_xCell_by_molecular_subtype.pdf`: heatmap of average immune scores per cell type per molecular subtype per histology.

### Running the analysis

```sh
bash run-immune-deconv.sh
```



