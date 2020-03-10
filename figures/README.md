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
|`na_color`|![f1f1f1](https://placehold.it/150x40/f1f1f1/FFFFFF?text=f1f1f1)| Single hex value|Values throughout the project that for various reasons are non-applicable|`color_key[is.na(color_key)] <- na_color`|
|`histologies_color_key`| Adenoma![0001FF](https://placehold.it/150x40/0001FF/FFFFFF?text=0001FF) <br> ATRT![001DCC](https://placehold.it/150x40/001DCC/FFFFFF?text=001DCC) <br> Central neurocytoma![006AFF](https://placehold.it/150x40/006AFF/FFFFFF?text=006AFF) <br> Chondrosarcoma![0071CC](https://placehold.it/150x40/0071CC/FFFFFF?text=0071CC) <br> Chordoma![00C5CC](https://placehold.it/150x40/00C5CC/FFFFFF?text=00C5CC) <br> Choroid plexus tumor![00CC2A](https://placehold.it/150x40/00CC2A/FFFFFF?text=00CC2A) <br> CNS EFT-CIC![00CC7E](https://placehold.it/150x40/00CC7E/FFFFFF?text=00CC7E) <br> CNS lymphoma![00D4FF](https://placehold.it/150x40/00D4FF/FFFFFF?text=00D4FF) <br> CNS neuroblastoma![00FF58](https://placehold.it/150x40/00FF58/FFFFFF?text=00FF58) <br> CNS Rhabdomyosarcoma![00FFC1](https://placehold.it/150x40/00FFC1/FFFFFF?text=00FFC1) <br> CNS sarcoma![12FF00](https://placehold.it/150x40/12FF00/FFFFFF?text=12FF00) <br> Craniopharyngioma![2ACC00](https://placehold.it/150x40/2ACC00/FFFFFF?text=2ACC00) <br> DNET![3800CC](https://placehold.it/150x40/3800CC/FFFFFF?text=3800CC) <br> Dysplasia![6691FF](https://placehold.it/150x40/6691FF/FFFFFF?text=6691FF) <br> Embryonal Tumor![66D0FF](https://placehold.it/150x40/66D0FF/FFFFFF?text=66D0FF) <br> Ependymoma![66FF70](https://placehold.it/150x40/66FF70/FFFFFF?text=66FF70) <br> ETMR![66FFB0](https://placehold.it/150x40/66FFB0/FFFFFF?text=66FFB0) <br> Ganglioglioma![66FFEF](https://placehold.it/150x40/66FFEF/FFFFFF?text=66FFEF) <br> Germinoma![6900FF](https://placehold.it/150x40/6900FF/FFFFFF?text=6900FF) <br> Glial-neuronal tumor NOS![7B66FF](https://placehold.it/150x40/7B66FF/FFFFFF?text=7B66FF) <br> Gliosis![7BFF00](https://placehold.it/150x40/7BFF00/FFFFFF?text=7BFF00) <br> Hemangioblastoma![7FCC00](https://placehold.it/150x40/7FCC00/FFFFFF?text=7FCC00) <br> Hemangioma![8C00CC](https://placehold.it/150x40/8C00CC/FFFFFF?text=8C00CC) <br> HGAT![9BFF66](https://placehold.it/150x40/9BFF66/FFFFFF?text=9BFF66) <br> Langerhans Cell histiocytosis![BA66FF](https://placehold.it/150x40/BA66FF/FFFFFF?text=BA66FF) <br> LGAT![CC00B8](https://placehold.it/150x40/CC00B8/FFFFFF?text=CC00B8) <br> LGMT![CC1C00](https://placehold.it/150x40/CC1C00/FFFFFF?text=CC1C00) <br> Medulloblastoma![CC7000](https://placehold.it/150x40/CC7000/FFFFFF?text=CC7000) <br> Meningioma![CCC500](https://placehold.it/150x40/CCC500/FFFFFF?text=CCC500) <br> MPNST![D200FF](https://placehold.it/150x40/D200FF/FFFFFF?text=D200FF) <br> Neurofibroma![DAFF66](https://placehold.it/150x40/DAFF66/FFFFFF?text=DAFF66) <br> none![E5FF00](https://placehold.it/150x40/E5FF00/FFFFFF?text=E5FF00) <br> Oligodendroglioma![f1f1f1](https://placehold.it/150x40/f1f1f1/FFFFFF?text=f1f1f1) <br> Other![F966FF](https://placehold.it/150x40/F966FF/FFFFFF?text=F966FF) <br> Pineoblastoma![FFA566](https://placehold.it/150x40/FFA566/FFFFFF?text=FFA566) <br> Schwannoma![FFB000](https://placehold.it/150x40/FFB000/FFFFFF?text=FFB000) <br> Teratoma![FFE566](https://placehold.it/150x40/FFE566/FFFFFF?text=FFE566) <br>|a named vector of the hex values that were assigned to each `short_histology` group. These are the same values as are seen in `palettes/histology_color_by_sample.tsv` but in a named vector format per histology format rather than as a per sample table|For color-coding by `short_histology` when its more convenient to provide a named vector| |
|`gradient_col_palette`| ![f7fcf5](https://placehold.it/150x40/f7fcf5/FFFFFF?text=f7fcf5) <br> ![e5f5e0](https://placehold.it/150x40/e5f5e0/FFFFFF?text=e5f5e0) <br> ![c7e9c0](https://placehold.it/150x40/c7e9c0/FFFFFF?text=c7e9c0) <br> ![a1d99b](https://placehold.it/150x40/a1d99b/FFFFFF?text=a1d99b) <br> ![74c476](https://placehold.it/150x40/74c476/FFFFFF?text=74c476) <br> ![41ab5d](https://placehold.it/150x40/41ab5d/FFFFFF?text=41ab5d) <br> ![238b45](https://placehold.it/150x40/238b45/FFFFFF?text=238b45) <br> ![006d2c](https://placehold.it/150x40/006d2c/FFFFFF?text=006d2c) <br> ![00441b](https://placehold.it/150x40/00441b/FFFFFF?text=00441b) <br>|Vector of 9 hex codes|For numeric data being plotted e.g. tumor mutation burden|
|`divergent_col_palette`|![67001f](https://placehold.it/150x40/67001f/FFFFFF?text=67001f) <br> ![b2182b](https://placehold.it/150x40/b2182b/FFFFFF?text=b2182b) <br> ![d6604d](https://placehold.it/150x40/d6604d/FFFFFF?text=d6604d) <br> ![f4a582](https://placehold.it/150x40/f4a582/FFFFFF?text=f4a582) <br> ![fddbc7](https://placehold.it/150x40/fddbc7/FFFFFF?text=fddbc7) <br> ![f7f7f7](https://placehold.it/150x40/f7f7f7/FFFFFF?text=f7f7f7) <br> ![d1e5f0](https://placehold.it/150x40/d1e5f0/FFFFFF?text=d1e5f0) <br> ![92c5de](https://placehold.it/150x40/92c5de/FFFFFF?text=92c5de) <br> ![4393c3](https://placehold.it/150x40/4393c3/FFFFFF?text=4393c3) <br> ![2166ac](https://placehold.it/150x40/2166ac/FFFFFF?text=2166ac) <br> ![053061](https://placehold.it/150x40/053061/FFFFFF?text=053061) <br>|Vector of 9 hex codes|For data has that is bidirectional e.g. Amplification/Deletion values like `seg.mean`|  |
|`binary_col_palette` |![67001f](https://placehold.it/150x40/67001f/FFFFFF?text=67001f) <br> ![053061](https://placehold.it/150x40/053061/FFFFFF?text=053061) <br>|A vector of two hex codes|For when we have two status values e.g. |    |  |
