## Expression Plots

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Purpose

1. Create `tumor vs normal` and `tumors only` expression plots, with `TPM` on the y-axis, as follows:

* `cohort + cancer_group tumor only plots`: For each cohort + cancer_group, plot all tumors.
* `cancer_group tumor only plots`: For each cancer_group, plot all tumors.
* `cohort + cancer_group tumor vs normal plots`: For each cancer_group within a cohort, plot cohort + cancer_group vs GTEx subgroups.
* `cancer_group tumor vs normal plots`: For each cancer_group, plot cancer_group vs GTEx subgroups.

2. Export two tables (.tsv) `cancer_group_level.tsv` and `cohort_cancer_group_level.tsv` with the following fields: `gene`, `ENSG_id`, `cohort`, `cancer_group`, `x_labels`, `mean`, `median`, `sd`, `efo_code`, `mondo_code`, `uberon_code`, `plot_api`.

3. Export metadata files to use in JSON files with the following fields: `gene` (gene name), `plot_type` (either `tumors_only` or `tumor_vs_normal`), `cohort` (`_` separated cohort name), `cancer_group` (only applicable for tumor vs normal plots and `NA` for tumors only plot), `analysis_type` (either `cohort_cancer_group_level` or `cancer_group_level`), `plot_fname` (plot filename) and `table_fname` (table filename).

### Methods 

1. Filter `cancer_subgroup` or `gtex_subgroup` with number of samples < 5 
2. Currently only generating outputs for `GPC2` and `MYCN` for the initial review

### Module structure

```
├── 01-tumor-gtex-plots.R # script to call all functions and run analysis
├── README.md 
├── plots # png files of expression boxplots (possible examples)
│   ├── GPC2_GMKF_Neuroblastoma_vs_GTEx_cohort_cancer_group_level.png
│   ├── GPC2_PBTA_Medulloblastoma_vs_GTEx_cohort_cancer_group_level.png
│   ├── GPC2_Ependymoma_vs_GTEx_cancer_group_level.png
│   ├── GPC2_PBTA_GMKF_pan_cancer_cancer_group_level.png
│   ├── GPC2_PBTA_GMKF_pan_cancer_cohort_cancer_group_level.png
│   └── ...
├── results
│   ├── cancer_group_level.tsv 
│   ├── cohort_cancer_group_level.tsv
│   └── metadata.tsv # mapping of filename and metadata for all comparisons 
├── run-tumor-gtex-plots.sh # full analysis script
├── run-tumor-only-plots.sh # script for tumor only plots
├── run-tumor-vs-normal-plots.sh # script for tumor vs gtex plots
└── util
    ├── pubTheme.R # publication quality ggplot2 theme
    ├── tumor_plot.R # function for tumor only plot
    └── tumor_vs_normal_plot.R # function for tumor vs normal plot
```

### Analysis scripts

#### 01-tumor-gtex-plots.R

##### Input parameters:

```
Rscript 01-tumor-gtex-plots.R --help

Options:
	--expr_mat=EXPR_MAT
		collapsed TPM expression data: HUGO gene symbol x Sample identifiers (.rds)

	--hist_file=HIST_FILE
		histologies file (.tsv)

	--cohort_list=COHORT_LIST
		comma separated list of cohorts

	--tumor_vs_normal=TUMOR_VS_NORMAL
		TRUE or FALSE

	--analysis_type=ANALYSIS_TYPE
		cohort_cancer_group_level or cancer_group_level

	--plot_width=PLOT_WIDTH
		width in pixels

	--plot_height=PLOT_HEIGHT
		height in pixels

	--mapping_file=MAPPING_FILE
		filename for writing out file names and other info
```

##### Inputs from data download:

Using the current version:

Expession matrix: `gene-expression-rsem-tpm-collapsed.rds` 
Histologies file: `histologies.tsv`

##### Outputs: 

1. `plots/*.png`: expression boxplots with `TPM` as y-axis

Tumors only filename format: 
`{gene_name}_{underscore_separated_cohort_name}_{analysis_type}.png` 

Tumors vs GTEx filename format: 
```
# cohort + cancer_group level
{gene_name}_{cohort_name}_{cancer_group_name}_vs_GTEx_{analysis_type}.png

# cancer_group level
{gene_name}_{cancer_group_name}_vs_GTEx_{analysis_type}.png
``` 

2. `results/cancer_group_level.tsv` and `results/cohort_cancer_group_level.tsv`: corresponding data for expression boxplots

3. `results/metadata.tsv`: mapping of filename and metadata for all comparisons 

### Running the analysis

```
# tumor only plots
bash run-tumor-only-plots.sh

# tumor vs normal plots
bash run-tumor-vs-normal-plots.sh

# running both scripts
bash run-tumor-gtex-plots.sh
```



