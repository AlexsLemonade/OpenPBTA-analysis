# GISTIC Comparision

**Module authors:** Chante Bethell ([@cbethell](https://github.com/cbethell)) and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

**Note: The files used to run the notebook in this directory were generated via [GISTIC](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/data-formats.md#gistic-output-file-formats) and stored in the `analyses/run-gistic/results` directory of this repository.
When re-running this module, you may want to regenerate the GISTIC output files using the most recent data release.**

## Usage

To run the R notebooks in this module from the command line, sequentially (assuming that you are in the top directory of this repository), use:

```
bash run-compare-gistic.sh
```

## Folder content

[`01-GISTIC-cohort-vs-histology-comparison.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/compare-gistic/01-GISTIC-cohort-vs-histology-comparison.nb.html) identifies, if any, disagreement between GISTIC results for the entire cohort versus the three largest histology groups (of the cohort) that we have GISTIC results for (LGAT, HGAT, and medulloblastoma).
The output of this notebook includes three multipanel plots, found in this `plots` directory, depicting the distribution of amplifications/deletions across chromosomes (using [G-scores](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-4-r41)) for the specified individual histology versus the entire PBTA cohort.

[`02-GISTIC-tidy-data-prep.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/compare-gistic/02-GISTIC-tidy-data-prep.nb.html) formats the [GISTIC files](https://www.genepattern.org/modules/docs/GISTIC_2.0) needed to compare GISTIC's copy number calls to the copy number calls we prepared in the [focal-cn-file-preparation](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/focal-cn-file-preparation/results) module.
**Note: This notebook will focus on preparing the GISTIC result files for the LGAT histology (the largest histology group in the cohort) and for the entire cohort. This decision was made for the purpose of developing this analysis in both a time efficient and informative manner.**
The output of this notebook includes tables to be consumed for both a gene-level comparison and cytoband-level comparison of copy number calls analysis notebooks.

Gene-level copy number call comparison tables are saved to  `cohort_gistic_gene_cn_status_table.tsv.gz` and `lgat_gistic_gene_cn_status_table.tsv.gz` and contain the following fields: 

| `gene_symbol` | `Kids_First_Biospecimen_ID` | `status` | `detection_peak` |
| ------------- | --------------------------- | ---------| ---------------- |

Cytoband-level copy number call comparison tables are saved to  `cohort_gistic_cytoband_cn_status_table.tsv.gz` and `lgat_gistic_cytoband_cn_status_table.tsv.gz` and contain the following fields: 

| `cytoband` | `Kids_First_Biospecimen_ID` | `status` |
| ---------- | --------------------------- | -------- |

## Folder Structure

```
├── 01-GISTIC-cohort-vs-histology-comparison.Rmd
├── 01-GISTIC-cohort-vs-histology-comparison.nb.html
├── 02-GISTIC-tidy-data-prep.Rmd
├── 02-GISTIC-tidy-data-prep.nb.html
├── README.md
├── plots
│   ├── hgat_gistic_scores_multipanel_plot.png
│   ├── lgat_gistic_scores_multipanel_plot.png
│   └── medulloblastoma_gistic_scores_multipanel_plot.png
├── results
│   ├── cohort_gistic_gene_cn_status_table.tsv.gz
│   ├── cohort_gistic_cytoband_cn_status_table.tsv.gz
│   ├── lgat_gistic_gene_cn_status_table.tsv.gz
│   └── lgat_gistic_cytoband_cn_status_table.tsv.gz
├── run-compare-gistic.sh
└── util
    └── GISTIC-comparison-functions.R

```
