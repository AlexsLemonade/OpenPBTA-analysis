## Publication Figures How to



### Running the figure creation script

To refresh all figures, call the main figures script:
```
bash scripts/run-figures.sh
```
This assumes your current directory is the top of this repository.

### Summary for each figure:

| Figure | Script | Run Requirements | Linked Analysis Modules |
|--------|--------|------------------|-------------------------|
|Figure 2 | [`scripts/fig2-mutational-landscape.R`](./scripts/fig2-mutational-landscape.R) | ~128MB of RAM are needed due to the run_caller_consensus_analysis-pbta.sh handling of large MAF files|[`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers)) [`mutational-signatures`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures)) |  
