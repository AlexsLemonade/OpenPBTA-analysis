# Mutation interaction plot (co-occurence and mutual exclusivity) generation

The scripts in this directory create a plot to display co-occurence and mutual exclusivity of mutations across tumors.

Currently this is done across all tumor types, with all available individuals for which whole genome sequences are available for at least one tumor.

Importantly, only a single sequencing sample from each individual is used. 
The analyses and plots created include information from the top 50 most mutated genes, with genes that are commonly mutated removed.
The commonly mutated genes are derived from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/, specifically the top 50 genes from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5706417/bin/12920_2017_309_MOESM3_ESM.txt

The main script creates plots for the full data set, as well as for groups of specific tumor types.
In the case of the full data set, a bar plot is also produced that summarizes which tumor types are mutated for each of the most common genes, as well as a publication-ready figure that combines the co-occurence plot and the bar plot.



### Example plot

![Co-occurence Plot](plots/consensus_top50.png)

## Usage

To run a basic analysis and create a co-occurence plot:
```
bash analyses/interaction-plots/01-create-interaction-plots.sh
```

# Result Ns used in manuscript
02-result-ns-for-manuscript.Rmd is a notebook which calculates the relevant Ns and ratios used in the manuscript.



