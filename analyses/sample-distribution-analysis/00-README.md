# OpenPBTA sample-distribution-analysis

## Folder content
This folder contains scripts tasked to analyze the distribution of samples across cancer types and histologies in the PBTA dataset.  
`00-install-packages.R` is a script written to install required packages and dependencies. This script is sourced in the other scripts.   
`01-pbta-analysis-5` is a script written to perform the initial analyses on the PBTA dataset (as discussed in [issue #5](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/5)). This script produces TSV files containing the counts and percentages of samples across each unique cancer type, and the tumor types sorted in order of descending expression at each unique brain location. 
`02-multilayer-pie.R` is a script written to produce an interactive treemap and multilayer pie chart representing the distribution of samples across broad histologies, short histologies, and tumor types.  
`03-histology-app.R` is a shiny script written to produce a Shiny web application displaying the multilayer pie chart.
`04-pbta-analysis-5.Rmd` is a R Markdown script written to summarize the results and plots produced within this folder. 
`04-pbta-analysis-5.nb.html` is the HTML output of the R Markdown mentioned above. 

## Folder structure 

The structure of this folder is as follows:

```
sample-distribution-analysis
├── 00-README.md
├── 00-install-packages.R
├── 01-pbta-analysis-5.R
├── 02-multilayer-pie.R
├── 03-histology-app.R
├── 04-pbta-analysis-5.nb.html
├── 04-pbta-analysis-5.Rmd
├── plots
│   ├── distribution_across_cancer_types.pdf
│   ├── histology-pie.html
│   └── histology-treemap.html
├── results
│   ├── disease_expression.tsv
│   ├── primary_sites.tsv
│   └── sunburst_plot_df.tsv

```