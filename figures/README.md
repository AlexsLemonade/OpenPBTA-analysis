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


## Color Palette Usage

This project has a unified color palette.
There are four sets of hex color keys to be used for all final figures.
They are created by running `scripts/color-palettes.R` which creates two files
in the `palettes` folder:
- `palettes/histology_color_by_sample.tsv`: A TSV file which contains hex codes for each
`Kids_First_Biospecimen_ID` based on their membership to `short_histology` from
the `data/pbta-histologies.tsv` file.
Biospecimens without a `short_histology` designation are assigned the `NA` color
for this project:
`#F1F1F1`

- `palettes/hex_color_palettes.rds`: contains a list of 5 vectors for color palettes:

| Palette Name | HEX codes included | Structure | Variable application | Example Usage |
|--------------|--------------------|-----------|----------------------|---------------|
|`na_color`|`#F1F1F1`| Single hex value|Values throughout the project that for various reasons are non-applicable|`color_key[is.na(color_key)] <- na_color`|
|`histologies_color_key`| See Table: `palettes/histology_color_by_sample.tsv` |a named vector of the hex values that were assigned to each `short_histology` group. These are the same values as are seen in `palettes/histology_color_by_sample.tsv` but in a named vector format per histology format rather than as a per sample table|For color-coding by `short_histology` when its more convenient to provide a named vector| |
|`gradient_col_palette`|`#f7fcf5` `#e5f5e0` `#c7e9c0` `#a1d99b` `#74c476` `#41ab5d` `#238b45` `#006d2c` `#00441b`|Vector of 9 hex codes|For numeric data being plotted e.g. tumor mutation burden|
|`divergent_col_palette`|`#67001f` `#b2182b` `#d6604d` `#f4a582` `#fddbc7` `#f7f7f7` `#d1e5f0` `#92c5de` `#4393c3` `#2166ac` `#053061`|Vector of 9 hex codes|For data has that is bidirectional e.g. Amplification/Deletion values like `seg.mean`|  |
|`binary_col_palette` |`#67001f` `#053061`|A vector of two hex codes|For when we have two status values e.g. |    |  |
