## Chromothripsis Analysis (in progress)

**Module authors :**
C. Savonen for ALSF CCDL, Krutika Gaonkar, Jo Lynne Rokita ([@jharenza](https://github.com/jharenza), and Run Jin ([@runjin326](https://github.com/runjin326) for D3B

This module runs survival analysis for groups or subgroups of samples based on certain histology categories.

#### Usage
The module can be run with the following bash script 
```
bash run_survival.sh
```

#### Inputs 

For `survival-analysis_HGG_DMG.Rmd` - all files are from data download
- `pbta-histologies.tsv`
- `independent-specimens.wgswxs.primary.tsv`

For `survival-analysis_histologies.Rmd` - 
- Two files from data download
  - `pbta-histologies.tsv`
  - `independent-specimens.rnaseq.primary-plus-stranded.tsv`
- Two more files from analysis module specified in the path are required:
  - `tp53_nf1_score/results/tp53_altered_status.tsv`
  - `telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_FPKM.txt`
- An additional file from `figures` folder is used for colors used in the figure:
  - `palettes/broad_histology_cancer_group_palette.tsv`

#### Scripts

Three notebooks are included in this analysis module:

- `survival-analysis_template.Rmd`
  - This notebook shows an example for running basic survival analysis models
  
- `survival-analysis_HGG_DMG.Rmd`
  - This notebook runs survival anaysis on all solid tumor DNA samples with `Diffuse astrocytic and oligodendroglial tumor` as `broad_histology`
  - Survival analysis for each `molecular_subtype` in this cohort is performed
  
- `survival-analysis_histologies.Rmd`
  - This notebook runs survival analysis on 3 different dataframes: 
    - all solid primary CNS tumor with RNA samples
    - all HGAT solid primary CNS tumor with RNA samples 
    - all non-HGAT solid primary CNS tumor with RNA samples with HGAT histology 
    
  - All samples are grouped into different categories based on TP53 and telomerase scores - the following categories (called `phenotypes` in this notebook) are included with the following rules
    - `tp53_0.5`: this columns indicates whether a sample is `tp53_high` or `tp53_low` based on whehter its TP53 score is over or below 0.5
    - `tel_0.5`: this columns indicates whether a sample is `telomerase_high` or `telomerase_low` based on whehter its telomerase score is over or below 0.5
    - `tp53_strict`: this columns indicates whether a sample is `tp53_high`, `tp53_mid`, or `tp53_low` based on whether they are above 75 quantile (`tp53_high`), below 25 quantile (`tp53_low`), or in between (`tp53_mid`)
    - `tel_strict`: this columns indicates whether a sample is `telomerase_high`, `telomerase_mid`, or `telomerase_low` based on whether they are above 75 quantile (`telomerase_high`), below 25 quantile (`telomerase_low`), or in between (`telomerase_mid`)
    - `pheno_0.5`: this column combines info from `tp53_0.5` and `tel_0.5`
    - `pheno_strict`: this column combines info from `tp53_strict` and `tel_strict`
    - `pheno_extremes`: this column only retains categorizations from `pheno_strict` column when the sample does not fit into `mid` group (between 25 and 75 percentile) for neither TP53 nor telomase scores
  
  -  Survival analyses were ran on the 3 different dataframes specified above, separated by the groups indicated in each column described above

#### Output tables
The following tables were output from each notebook
- `survival-analysis_template.Rmd`
  - `results/cox_regression_tmb.tsv`
  - `results/logrank_gender.tsv`

- `survival-analysis_HGG_DMG.Rmd`
  - `results/survival_model_DMG_H3_K28_v_HGG_H3_wildtype.tsv`

- `survival-analysis_histologies.Rmd`
  - The models for each phenotype in all three different dataframes were output to `results/model` folder as `{data_frame}-{phenotype}-model.tsv` 
  - The statistics for each phenotype in all three different dataframes were output to `results/stats` folder as `adjusted_stats_for_{data_frame}_by_{phenotype}.txt` 


#### Plots
The following plots were output from each notebook
- `survival-analysis_template.Rmd`
  - `plots/survival_curve_gender.pdf`

- `survival-analysis_HGG_DMG.Rmd`
  - `plots/survival_curve_DMG_H3_K28_v_HGG_H3_wildtype.pdf`

- `survival-analysis_histologies.Rmd`
  - The survival plots for each phenotype in all three different dataframes were output to `plots` folder as
`{data_frame}-{phenotype}-KM.pdf`


