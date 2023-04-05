
# Data underlying OpenPBTA Manuscript Figures and Molecular Alterations

This directory contains CSV files that represent data contained in plots shown in the OpenPBTA manuscript.
It is intended to facilitate inspection of the underlying data shown in each figure and to explicitly capture which samples are included in figures (where applicable).
To **reproduce the figures**, we recommend using the code in the analysis repository: <https://github.com/AlexsLemonade/OpenPBTA-analysis>.
Please see the `figures/` directory documentation in the repository and the documentation for figure generation scripts (`figures/scripts/README.md`).


<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Figure Data](#figure-data)
  - [Main Figure Data](#main-figure-data)
    - [Figure 1](#figure-1)
    - [Figure 2](#figure-2)
    - [Figure 3](#figure-3)
    - [Figure 4](#figure-4)
    - [Figure 5](#figure-5)
  - [Supplemental Figure Data](#supplemental-figure-data)
    - [Figure S2](#figure-s2)
    - [Figure S3](#figure-s3)
    - [Figure S4](#figure-s4)
    - [Figure S5](#figure-s5)
    - [Figure S6](#figure-s6)
    - [Figure S7](#figure-s7)
- [Molecular Alterations](#molecular-alterations)
  - [File Description](#file-description)
  - [Alteration Notation](#alteration-notation)
  - [Caveats](#caveats)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Figure Data

Data for the figures below are included as part of this download.
The table for each figure contains the following:

| Column | Purpose |
|--------|---------|
| Figure Panel | Number and letter of the panel(s) |
| Filename | Name of the file that contains the relevant panel's data |
| Brief Description | An overview of what the panel includes; please see the manuscript [figure legends](https://alexslemonade.github.io/OpenPBTA-manuscript/#figure-titles-and-legends) and [supplemental figures](https://alexslemonade.github.io/OpenPBTA-manuscript/#supplemental-information-titles-and-legends) for more information |
| Notes on Tabular Data | Any notes to aid downloaders in interpreting the table; may be left blank |

Where applicable, all tables are ordered by either `sample_id` or `Kids_First_Biospecimen_ID`.


### Main Figure Data

#### Figure 1

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   1B     | `figure-1b-data.csv`  | Barplot showing the number of biospecimens per phase of therapy | |

#### Figure 2

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   2A     | `figure-2a-data.csv`   | Oncoprint summarizing low-grade glioma alterations across samples for top mutated genes           |
|   2B     | `figure-2b-data.csv`   | Oncoprint summarizing embryonal tumor alterations across samples for top mutated genes         |
|   2C     | `figure-2c-data.csv`   | Oncoprint summarizing high-grade glioma alterations across samples for top mutated genes         |
|   2D     | `figure-2d-data.csv`   |  Oncoprint summarizing other CNS tumor alterations across samples for top mutated genes         |
#### Figure 3

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   3A     | `figure-3a-data.csv`  | Barplot of occurrence and co-occurrence of nonsynonymous mutations for the 50 most commonly mutated genes across all tumor types         |
|   3B     | `figure-3b-data.csv`  | Heatmap of co-occurrence and mutual exclusivity of nonsynonymous mutations between genes           |
|   3C     | `figure-3c-data.csv`  | Scatterplot showing correlation of SV and CNV breaks            |
|   3D     | `figure-3d-data.csv`  | Barplot showing chromothripsis frequency for all cancer groups with N >= 3 tumors           |
|   3E     | `figure-3e-data.csv`  | Sina plots of RefSig signature weights across cancer groups         |
#### Figure 4

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   4A     | `figure-4a-data.csv` | ROC curve for TP53 classifier run on stranded RNA-Seq data     |
|   4B     | `figure-4b-data.csv` | Violin/strip plots of TP53 scores across TP53 alteration status           |
|   4C     | `figure-4c-data.csv` |  Violin/strip plots of TP53 expression [log(FPKM)] across TP53 alteration status            |
|   4D     | `figure-4d-data.csv` | Box/strip plots of TP53 scores and telomerase scores across cancer groups          |
|   4E     | `figure-4e-data.csv` |  Heatmap of RefSig signatures in patients with a hypermutator phenotype       |
|   4F     | `figure-4f-data.csv`  | Forest plot from survival analysis exploring TP53 and telomerase score effects          |
|   4G     | `figure-4g-data.csv`  | Forest plot from survival analysis exploring HGG molecular subtype effects          |
|   4H     | `figure-4h-data.csv` | Kaplan-Meier cure of HGG tumors by molecular subtype           |

#### Figure 5

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   5A     | `figure-5a-data.csv`  | UMAP of tumors colored by broad histology           |
|   5B     | `figure-5b-data.csv`  | Heatmap of GSVA scores across cancer groups           |
|   5C     | `figure-5c-data.csv`  | Box/strip plots of quanTIseq scores across cancer groups           |
|   5D     | `figure-5d-data.csv`  | Forest plot from survival analysis exploring CD274 expression and immune cell proportion effects           |


### Supplemental Figure Data

#### Figure S2

| Figure Panel(s) | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   S2A-C    |  `figure-S2a-S2b-S2c-data.csv`   | Correlations, distributions, and UpSet plot of mutation VAFs from variant callers for PBTA data, respectively   | This compressed file contains data for all three panels S2A, S2B, and S2C. Caution: Uncompressed, this file is ~3.9 GB.|
|   S3D-F     |   `figure-S2d-S2e-S2f-data.csv`    |  Correlations, distributions, and UpSet of mutation VAFs from variant callers for TCGA data, respectively       | This compressed file contains data for all three panels S2D, S2E, and S2F. Uncompressed, this file is 13MB.            |
|   S2G     |   `figure-S2g-data.csv`          |  Distributions of VAFs from Lancet calls on PBTA data for tumors with WGS and WXS   |
|   S2H     |   `figure-S2h-data.csv`         |  Cumulative distribution TMB plots for PBTA cancer groups          |
|   S2I     |   `figure-S2i-data.csv`         |  Cumulative distribution PBTA plots for TCGA histologies          |

#### Figure S3

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   S3A     |    `figure-S3a-data.csv`         |   Violin plots of tumor purity by cancer group          |
|   S3B     |    `figure-S3b-data.csv`          |  Oncoprint summarizing rare CNS tumor alterations across samples for top mutated genes    |
|   S3C     |    `figure-S3c-data.csv.gz`       | Genome-wide plot of CNV alterations for each sample, grouped by broad histology  |   This file is compressed; uncompressed, this file is ~15 MB.      |
|   S3D     |    `figure-S3d-data.csv`          | Box plots of SNV breaks across number of chromothripsis regions            |
|   S3E     |    `figure-S3e-data.csv`          | Box plots of CNV breaks across number of chromothripsis regions            |


#### Figure S4

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   S4A     |   `figure-S4a-data.csv`           |  Sample-specific signature weights for 8 adult CNS-specific RefSig mutational signatures           |
|   S4B     |   `figure-S4b-data.csv`           |   RefSig Signature 1 weights across cancer groups by phrase of therapy          |


#### Figure S5

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   S5A     |  `figure-S5a-data.csv`            |  ROC curve for _TP53_ classifier run on poly-A RNA-Seq samples           |
|   S5B     |  `figure-S5b-data.csv`            |  Correlation of _TP53_ classifier score with _TERT_ expression (FPKM)           |
|   S5C     |  `figure-S5c-data.csv`            |  Correlation of _TP53_ classifier score with _TERC_ expression (FPKM)           |           |



#### Figure S6

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   S6A     |   `figure-S6a-data.csv`          |  UMAP of stranded RNA-Seq samples, with medulloblastoma molecular subtypes emphasized           |
|   S6B     |   `figure-S6b-data.csv`          |  UMAP of stranded RNA-Seq samples, with ependymoma molecular subtypes emphasized           |
|   S6C     |   `figure-S6c-data.csv`          |   UMAP of stranded RNA-Seq samples, with low-grade glioma molecular subtypes emphasized          |
|   S6D     |   `figure-S6d-data.csv`          |  UMAP of stranded RNA-Seq samples, with high-grade glioma molecular subtypes emphasized            |
|   S6E     |   `figure-S6e-data.csv`          |   quanTIseq-estimated immune cell fractions in histologies with >=3 molecular subtypes         |
|   S6F     |   `figure-S6f-data.csv`          |  Ratio of CD8+/CD4+ T-cell fractions in histologies with >=3 molecular subtypes             |


#### Figure S7

| Figure Panel | Filename | Brief Description | Notes on Tabular Data |
|--------|-------------|--------------------|-------------------------|
|   S7A     |   `figure-7a-data.csv`           |  Counts of RNA-Seq samples across cancer groups colored by library preparation method           |
|   S7B     |   `figure-7b-data.csv`           |   UMAP of all RNA-Seq samples, colored by cancer group, with shapes indicating library preparation method           |
|   S7C     |   `figure-7c-data.csv`           |  UMAP of stranded RNA-Seq samples, colored by cancer group, with shapes indicating sequencing center            |
|   S7D     |   `figure-7d-data.csv`           |   ROC curve for TP53 classifier run on stranded RNA-Seq data, considering only tumors with high tumor purity     | High tumor purity is defined as having a tumor purity >= median of the tumor's respective cancer group. <br> This figure corresponds to 4A in the manuscript.      |
|   S7E     |   `figure-7e-data.csv`           |  Violin/strip plots of TP53 scores across TP53 alteration status, considering only tumors with high tumor purity            | High tumor purity is defined as having a tumor purity >= median of the tumor's respective cancer group. <br> This figure corresponds to 4B in the manuscript.           |
|   S7F     |   `figure-7f-data.csv`           | Violin/strip plots of TP53 expression [log(FPKM)] across TP53 alteration status, considering only tumors with high tumor purity             |   High tumor purity is defined as having a tumor purity >= median of the tumor's respective cancer group. <br> This figure corresponds to 4C in the manuscript.         |
|   S7G     |   `figure-7g-data.csv`           |  UMAP of tumors colored by broad histology, considering only tumors with high tumor purity            |  High tumor purity is defined as having a tumor purity >= median of the tumor's respective cancer group. <br> This figure corresponds to 5A in the manuscript.          |
|   S7H     |   `figure-7h-data.csv`           |  Box/strip plots of quanTIseq scores across cancer groups, considering only tumors with high tumor purity           |  High tumor purity is defined as having a tumor purity >= median of the tumor's respective cancer group. <br> This figure corresponds to 5C in the manuscript.           |
|   S7I     |   `figure-7i-data.csv`           |  Box/strip plots of TP53 scores and telomerase scores across cancer groups, considering only tumors with high tumor purity       | High tumor purity is defined as having a tumor purity >= median of the tumor's respective cancer group. <br> This figure corresponds to 4D in the manuscript.          |


## Molecular Alterations

Please see the `tables/tabulate-molecular-alterations.R` script in the `OpenPBTA-analysis` repository for more information about how this table was generated.

### File Description

`openpbta-molecular-alterations.csv` contains information about the alterations (SNV, CNV, and fusions) a sample has in any gene included in the oncoprints in the manuscript (i.e., Figure 2 and Figure S3B).

In the case of multiple alterations affecting the same gene, individual alterations are separated by semi-colons.
We use the value `None` when no gene alterations are detected in the consensus data.

Each row corresponds to a `sample_id`-`composition` pair in the OpenPBTA data, where a given `sample_id` can be associated with multiple `composition` values (here, solid tissue or derived cell line).

We include `Kids_First_Biospecimen_DNA` and `Kids_First_Biospecimen_RNA` columns with the biospecimen IDs for DNA (WGS, WXS, Targeted) and RNA (RNA-seq) assays, respectively.

### Alteration Notation

* For SNV data, we prioritize values from the consensus MAF file (`pbta-snv-consensus-mutation.maf.tsv.gz`) for inclusion in the following order: `HGVSp_Short`, `HGVSc`, and `Variant_Type`.
When no change in the protein is noted in the `HGVSp_Short` value, we use the nucleotide change.
* CNV alterations use the following notation from the `consensus_seg_annotated_*` files included in the data download: `<cytoband>-<status>-<copy_number>`.
* We report the `FusionName` field from `pbta-fusion-putative-oncogenic.tsv` for fusions.
When reciprocal fusions are detected, we only report one – whichever comes first when partner genes are sorted alphabetically.

### Caveats

* Some `sample_id`-`composition` pairs map to multiple biospecimens of the same type (i.e., DNA and RNA).
When this occurs, the individual biospecimen identifiers are separated with semi-colons in the relevant column, and `multiple_assays_within_type` is marked as `TRUE`.
The alterations for a gene are the union of alterations detected in biospecimens (i.e., they are detected in at least one biospecimen).
* If multiple `cancer_type` or `broad_histology` values are associated with the same `sample_id`-`composition` pair, these values are semi-colon separated.
* `germline_sex_estimate` is only available for `sample_id`-`composition` pairs that were assayed with WGS.
* In two instances, multiple RNA biospecimens map to the sample `sample_id`-`composition` pair, and only one RNA biospecimen was assigned the `germline_sex_estimate` derived from the WGS data.
We only included the non-missing values in these rows.
The two `sample_id` instances impacted are `7316-161` and `7316-85`.
