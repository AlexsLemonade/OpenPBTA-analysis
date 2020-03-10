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
All figures print out to the `figures/pngs` folder and will be linked to the
accompanying manuscript repository `OpenPBTA-manuscript`.

## Summary for each figure:

Each figure has its own script stored in the `figures/scripts`.
All are called by the main bash script, `figures/run-figures.sh`.
However, some are have higher memory requirements than others.

| Figure | Script | Run Requirements | Linked Analysis Modules | Files consumed |
|--------|--------|------------------|-------------------------|----------------|
|Figure 2 | [scripts/fig2-mutational-landscape.R`](./scripts/fig2-mutational-landscape.R) | ~128MB of RAM are needed due to the run_caller_consensus_analysis-pbta.sh handling of large MAF files|[`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers)) [`mutational-signatures`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures)) |  `pbta-snv-lancet.vep.maf.gz` <br> `pbta-snv-mutect2.vep.maf.gz` <br> `pbta-snv-strelka2.vep.maf.gz` <br> `pbta-snv-vardict.vep.maf.gz` <br> `tcga-snv-lancet.vep.maf.gz` <br> `tcga-snv-mutect2.vep.maf.gz` <br> `tcga-snv-strelka2.vep.maf.gz` |
