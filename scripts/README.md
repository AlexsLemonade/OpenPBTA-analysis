## Running the figure creation script

Steps to refreshing all publication-ready figures, these steps assumes your
current directory is the top of this repository. 

1) Have the current dataset. 
See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#how-to-obtain-openpbta-data)
about obtaining the current data release. 

2) Set up the Docker container
See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#docker-image) 
about setting up this repository's Docker container.

3) Call the main figures bash script
This script will run all the steps needed to take the
```
bash scripts/run-figures.sh
```
All figures print out to the `./figures` folder and will be linked to the 
accompanying manuscript repository `OpenPBTA-manuscript`.
 
4) Refresh manubot? 


### Summary for each figure:

| Figure | Script | Run Requirements | Linked Analysis Modules |
|--------|--------|------------------|-------------------------|
|Figure 2 | [`scripts/fig2-mutational-landscape.R`](./scripts/fig2-mutational-landscape.R) | ~128MB of RAM are needed due to the run_caller_consensus_analysis-pbta.sh handling of large MAF files|[`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers)) [`mutational-signatures`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures)) |  
