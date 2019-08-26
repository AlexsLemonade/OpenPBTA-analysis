# OpenPBTA sample-distribution-analysis

## Usage

To run the scripts in this module from the command line, use:

```
bash "analyses/sample-distribution-analysis/run-sample-distribution.sh"

```
## Folder content
This folder contains scripts tasked to analyze the distribution of samples across cancer types and histologies in the PBTA dataset.

`01-filter-across-types.R` is a script written to perform the initial analyses on the PBTA dataset (as discussed in [issue #5](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/5)).  
This script produces TSV files containing the counts and percentages of samples across each unique cancer type, and the types sorted in order of descending expression at each unique brain location.   

`02-multilayer-pie.R` is a script written to produce an interactive treemap and multilayer pie chart representing the distribution of samples across broad histologies, short histologies, and tumor types.  

 The interactive treemap can be found [here](alexslemonade.github.io/OpenPBTA-analysis/histology-treemap.html).  
 The interactive, multilayer pie chart can be found [here](alexslemonade.github.io/OpenPBTA-analysis/histology-pie.html).

## Folder structure 

The structure of this folder is as follows:

```
sample-distribution-analysis
├── README.md
├── 01-filter-across-types.R
├── 02-multilayer-plots.R
├── plots
│   ├── distribution_across_cancer_types.pdf
│   ├── histology-pie.html
│   └── histology-treemap.html
├── results
│   ├── disease_expression.tsv
│   ├── primary_sites.tsv
│   └── sunburst_plot_df.tsv

```