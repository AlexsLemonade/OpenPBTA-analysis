# OpenPBTA sample-distribution-analysis

## Usage

To run all of the Rscripts in this module from the command line sequentially, use:

```
bash run-sample-distribution.sh
```

`run-sample-distribution.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

## Folder content

This folder contains scripts tasked to analyze the distribution of samples across cancer types and histologies in the PBTA dataset.

`01-filter-across-types.R` is a script written to perform the initial analyses on the PBTA dataset (as discussed in [issue #5](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/5)).  
This script produces TSV files containing the counts and percentages of samples across each unique cancer type, and the types sorted in order of descending expression at each unique brain location.   

`02-multilayer-pie.R` is a script written to produce an interactive treemap and multilayer pie chart representing the distribution of samples across broad histologies, short histologies, and tumor types.  
 
 The interactive treemap can be found [here](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/sample-distribution-analysis/plots/histology-treemap.html).  
 The interactive, multilayer pie chart can be found [here](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/plots/sample-distribution-analysis/histology-pie.html).

The `03-tumor-descriptor-and-assay-count` notebook contains a series of tables that count the number of each assay type and example primary vs. recurrence broken down by histology. 
View the notebook [here](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/sample-distribution-analysis/03-tumor-descriptor-and-assay-count.nb.html).



## Folder structure 

The structure of this folder is as follows:

```
├── 01-filter-across-types.R
├── 02-multilayer-plots.R
├── 03-tumor-descriptor-and-assay-count.Rmd
├── 03-tumor-descriptor-and-assay-count.nb.html
├── README.md
├── plots
│   ├── distribution_across_cancer_types.pdf
│   ├── histology-pie.html
│   └── histology-treemap.html
├── results
│   ├── disease_expression.tsv
│   ├── plots_df.tsv
│   └── primary_sites_counts.tsv
└── run-sample-distribution.sh
```
