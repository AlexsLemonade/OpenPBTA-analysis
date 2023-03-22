This module creates tables and supplementary tables for the manuscript, as well as for upload to Zenodo.
All tables represent results from the most current OpenPBTA release.

Before creating tables, you must ensure that all manuscript analysis and figure generation scripts have been run.
For documentation on this, please see [the `scripts/` directory `README`](../scripts/README.md) and the [`figures/` directory `README`](../figures/README.md).


To run this module, run:

```
bash run-tables.sh
```

This will execute, in order,
- `write-manuscript-tables.Rmd`, which will export tables to one of `manuscript-tables/` or `manuscript-tables/supp/`, for main text and supplementary tables respectively.
Note that this consumes `input/` files to create certain tables.
- `util/BS_F0GNWEJJ_genomic_investigation.Rmd` which performs an investigation into a particular hypermutator sample of interest present in Table 2.
- `copy-zenodo-tables.sh`, which will copy any figure data tables exported within `analysis/` modules into `zenodo-upload/`.

## Manuscript tables

The following tables are exported from `write-manuscript-tables.Rmd`:

Table 1: Molecular subtypes determined for this project
- This table shows the type and number of different molecular subtypes profiled in the project.

Table 2: Patients with hypermutant tumors
- Listed are patients with at least one hypermutant or ultra-hypermutant tumor (inclusive of derived cell lines and all phases of therapy).
Coding region TMB, phase of therapy, therapeutic interventions, cancer predispositions, pathogenic germline variants, and molecular subtypes are included.


Table S1: Histologies table
- This table displays the histologies file in an excel format.
    - Sheet 1: README
    - Sheet 2: histologies file
    - Sheet 3: CNS region definitions

Table S2: DNA results table
- This excel shows key results from analyzing DNA sequencing.
    - Sheet 1: this table shows tumor mutation burden (TMB) for all samples - `TMB_coding` is calculated based on mutations in coding regions and `TMB_all` is calculated by all mutations a sample has
        - Input files `pbta-snv-consensus-mutation-tmb-all.tsv` and `pbta-snv-consensus-mutation-tmb-coding.tsv` are in data release.
    - Sheet 2: this table shows signature weights for CNS signature identified in [Degasperi et al, 2020](https://doi.org/10.1038/s43018-020-0027-5) using [`deconstructSigs`](https://doi.org/10.1186/s13059-016-0893-4) per sample
        - Input file `deconstructsigs_exposures_merged.tsv` is generated from module [mutational-signatures](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures)
    - Sheet 3: this table summarizes the number of chromothripsis events per sample
        - Input file `chromothripsis_summary_per_sample.txt` is generated from module [chromothripsis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromothripsis)

Table S3: RNA results table
- This excel shows key results from analyzing RNA sequencing.
    - Sheet 1: this table shows TP53 scores, cancer prediposition, SNV indel counts, CNV loss counts, SV counts, fusion counts and details about these events for each sample
        - Input file `tp53_altered_status.tsv` is from module [tp53_nf1_scores](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tp53_nf1_score)
    - Sheet 2: this tables shows the EXTEND scores for all samples - `NormEXTENDScores_counts` is EXTEND scores calculated based on expression counts matrix and `NormEXTENDScores_fpkm` is EXTEND scores calculated based on expression FPKM matrix
        - Input files are from module [telomerase-activity-prediction](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction)
    - Sheet 3: this tables shows the fractions of for various immune cell types in the tumor microenvironment (TME) for each RNA sample as calculated with [`quanTIseq`](https://doi.org/10.1186/s13073-019-0638-6).
        - Input files are from module [immune-deconv](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/immune-deconv)

Table S4: Overall survival models
- Table of overall survival analyses included in the manuscript.
    - Sheet 1: cox regression ~ `molecular_subtype` for HGG
    - Sheet 2: log-rank ~ `molecular_subtype` for HGG
    - Sheet 3: cox regression ~ `tp53 scores + telomerase scores + extent of tumor resection + LGG group + HGG group` across PBTA
    - Sheet 4: cox regression ~ `tp53 scores + telomerase scores + extent of tumor resection` for EPN
    - Sheet 5: cox regression ~ `tp53 scores + telomerase scores + extent of tumor resection` for DMG
    - Sheet 6: cox regression ~ `quantiseq cell type fractions + CD274 expression + extent of tumor resection` for MB

Table S5: Key Resources
- Table of all software and their respective versions used for the OpenPBTA project. Of note, this table contains all software in the OpenPBTA docker image utilized within the repository, but not all software was used for the final manuscript.
    - Sheet 1: r_packages
    - Sheet 2: python_libraries
    - Sheet 3: other_command_line_tools
    - Sheet 4: workflow_repository_tools


## Zenodo Upload

Files in this directory are named as: `figure-<figure number><figure panel>-data.csv` (or `.csv.gz` if compressed.)
For example, `figure-3a-data.csv` contains data for Figure 3A.
