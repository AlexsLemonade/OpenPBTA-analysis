## Running the figure generation script

Steps to refresh all publication-ready figures. 
All steps assume your current directory is the top of this repository.

#### 1. Obtain the current dataset.

See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#how-to-obtain-openpbta-data) to obtain the current data release.
We recommend [using the download script](https://github.com/AlexsLemonade/OpenPBTA-analysis#data-access-via-download-script) to obtain data because this will automatically create symlinks in `data/` to the latest files. 

#### 2. Set up an up-to-date project Docker container.

See [these instructions](https://github.com/AlexsLemonade/OpenPBTA-analysis#docker-image) for setting up the project Docker container.
Briefly, the latest version of the project Docker image, which is updated upon commit to `master`, can be obtained and run via:
```
docker pull ccdlopenpbta/open-pbta:latest
docker run \
  -e PASSWORD=<password> \
  -p 8787:8787 \
  -v $(pwd):/home/rstudio/kitematic \
  ccdlopenpbta/open-pbta:latest
```
You may choose to use [`docker exec`](https://docs.docker.com/engine/reference/commandline/exec/) to interact with the container from there or if you'd prefer the RStudio interface, you can navigate to `localhost:8787` and enter username `rstudio` and the password you set with the `run` command above.

#### 3. Run the bash script that generates the figures (`scripts/run-figures.sh`). 

This script runs **_all_** the intermediate steps needed to generate figures starting with the original data files.

```
bash scripts/run-figures.sh
```

Figures are saved to the `figures/pngs` folder and will be linked to the accompanying manuscript repository [`AlexsLemonade/OpenPBTA-manuscript`](https://github.com/AlexsLemonade/OpenPBTA-manuscript/).

## Summary for each figure

Each figure has its own script stored in the `figures/scripts`.
All are called by the main bash script `figures/run-figures.sh`.
However, we list information about the resources, intermediate steps, and [PBTA data files](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/data-formats.md#pbta-data-files) required for generating each figure below for convenience.

| Figure | Individual script | Notes on requirements | Linked analysis modules | PBTA data files consumed |
|--------|--------|------------------|-------------------------|-----------------------------|
| Figure 2 | [scripts/fig2-mutational-landscape.R`](./scripts/fig2-mutational-landscape.R) | ~128MB of RAM are needed due to the run_caller_consensus_analysis-pbta.sh handling of large MAF files|[`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) <br> [`mutational-signatures`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures) |  `pbta-snv-lancet.vep.maf.gz` <br> `pbta-snv-mutect2.vep.maf.gz` <br> `pbta-snv-strelka2.vep.maf.gz` <br> `pbta-snv-vardict.vep.maf.gz` <br> `tcga-snv-lancet.vep.maf.gz` <br> `tcga-snv-mutect2.vep.maf.gz` <br> `tcga-snv-strelka2.vep.maf.gz` |
