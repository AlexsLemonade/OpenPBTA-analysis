## Expression Plots

**Module authors:** Komal Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Purpose

1. Create `tumor-normal-gtex` and `pan-cancer` expression plots, with `TPM` on the y-axis, as follows:

* `cohort + cancer_group pan-cancer plots`: For each cohort + cancer_group, plot all tumors.
* `cancer_group pan-cancer plots`: For each cancer_group, plot all tumors.
* `cohort + cancer_group tumor-normal-gtex plots`: For each cancer_group within a cohort, plot cohort + cancer_group vs GTEx subgroups.
* `cancer_group tumor-normal-gtex plots`: For each cancer_group, plot cancer_group vs GTEx subgroups.

2. Export two tables (.tsv) `cancer_group_level.tsv` and `cohort_cancer_group_level.tsv` per comparison type i.e. `tumor_normal_gtex_plots_` and `pan_cancer_plots_` with the following fields: `Gene_symbol`, `Gene_Ensembl_ID`, `Dataset`, `Disease`, `x_labels`, `mean`, `median`, `sd`, `plot_api`. Add `MONDO`, `RMTL` and	`EFO` using `long-format-table-utils/annotator/annotator-cli.R`. Convert to JSONL and gzip.

3. Export metadata files with the following fields: `gene` (gene name), `plot_type` (either `pan_cancer` or `tumor_normal_gtex`), `cohort` (tumor cohort name for tumor-normal-gtex plots and `all_cohorts` for pan-cancer plots) for , `cancer_group` (only applicable for tumor-normal-gtex plots and `NA` for tumors only plot), `analysis_type` (either `cohort_cancer_group_level` or `cancer_group_level`), `plot_fname` (plot filename) and `table_fname` (table filename).

### Methods 

1. Filter `cancer_subgroup` or `gtex_subgroup` with number of samples < 5 

### Module structure

```
├── 01-tumor-gtex-plots.R # script to call all functions and run analysis
├── 02-chunkwise-annotate-and-zip.R # script to annotate files with MONDO, RMTL, EFO fields, convert to JSONL and gzip
├── README.md 
├── plots # png files of expression boxplots (possible examples)
│   ├── GPC2_GMKF_Neuroblastoma_vs_GTEx_cohort_cancer_group_level.png
│   ├── GPC2_PBTA_Medulloblastoma_vs_GTEx_cohort_cancer_group_level.png
│   ├── GPC2_Ependymoma_vs_GTEx_cancer_group_level.png
│   ├── GPC2_PBTA_GMKF_pan_cancer_cancer_group_level.png
│   ├── GPC2_PBTA_GMKF_pan_cancer_cohort_cancer_group_level.png
│   └── ...
├── results
│   ├── metadata.tsv # metadata file for all comparisons
│   ├── pan_cancer_plots_cancer_group_level.jsonl.gz
│   ├── pan_cancer_plots_cancer_group_level.tsv
│   ├── pan_cancer_plots_cohort_cancer_group_level.jsonl.gz
│   ├── pan_cancer_plots_cohort_cancer_group_level.tsv
│   ├── tumor_normal_gtex_plots_cancer_group_level.jsonl.gz
│   ├── tumor_normal_gtex_plots_cancer_group_level.tsv
│   ├── tumor_normal_gtex_plots_cohort_cancer_group_level.jsonl.gz
│   └── tumor_normal_gtex_plots_cohort_cancer_group_level.tsv 
├── run-pan-cancer-plots.sh # script for pan-cancer plots
├── run-tumor-gtex-plots.sh # full analysis script
├── run-tumor-normal-gtex-plots.sh # script for tumor vs gtex plots
├── run-annotator.sh # script to call 02-chunkwise-annotate-and-zip.R
└── util
    ├── pubTheme.R # publication quality ggplot2 theme
    ├── pan_cancer_plot.R # function for pan-cancer plot
    └── tumor_normal_gtex_plot.R # function for tumor-normal-gtex plot
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

	--map_file=MAP_FILE
		gene symbol-ensembl id mapping file

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

	--meta_file=META_FILE
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

2. `results/pan_cancer_plots_cancer_group_level.tsv` and `results/pan_cancer_plots_cohort_cancer_group_level.tsv`: corresponding data for pan-cancer expression boxplots

3. `results/tumor_normal_gtex_plots_cancer_group_level.tsv` and `results/tumor_normal_gtex_plots_cohort_cancer_group_level.tsv`: corresponding data for tumor-normal-gtex expression boxplots

4. `results/metadata.tsv`: metadata file for all comparisons 

#### 02-chunkwise-annotate-and-zip.R

##### Input parameters:

```
Rscript 02-chunkwise-annotate-and-zip.R --help

Options:
	--input_file=INPUT_FILE
		input file to annotate, convert to jsonl and gzip
```

##### Inputs:

This script takes the output `.tsv` files from `01-tumor-gtex-plots.R`, annotates using the 
`long-format-table-utils/annotator/annotator-cli.R` script, reads back annotated data (.tsv) into chunks of 1 million rows, converts to JSON and then JSONL, followed by combining all JSONL files and gzipping the final JSONL file. 

##### Outputs:

1. `results/pan_cancer_plots_cancer_group_level.jsonl.gz` and `results/pan_cancer_plots_cohort_cancer_group_level.jsonl.gz`: corresponding data for pan-cancer expression boxplots

2. `results/tumor_normal_gtex_plots_cancer_group_level.jsonl.gz` and `results/tumor_normal_gtex_plots_cohort_cancer_group_level.jsonl.gz`: corresponding data for tumor-normal-gtex expression boxplots

### Running the analysis

```
# pan-cancer plots
bash run-pan-cancer-plots.sh

# tumor-normal-gtex plots
bash run-tumor-normal-gtex-plots.sh

# annotator script
bash annotator.sh

# running both scripts
bash run-tumor-gtex-plots.sh
```

### Final output

The final plots as `.png` and tables as `.tsv` and `.jsonl.gz` files have been uploaded to `s3://kf-openaccess-us-east-1-prd-pbta/open-targets/analyses/tumor-gtex-plots/`



