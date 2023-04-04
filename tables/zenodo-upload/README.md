
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
|   1B     |  | Barplot showing the number of biospecimens per phase of therapy | |

#### Figure 2

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   2A     |   | Oncoprint summarizing low-grade glioma alterations across samples for top mutated genes           |
|   2B     |    | Oncoprint summarizing embryonal tumor alterations across samples for top mutated genes         |
|   2C     |    | Oncoprint summarizing high-grade glioma alterations across samples for top mutated genes         |
|   2D     |   |  Oncoprint summarizing other CNS tumor alterations across samples for top mutated genes         |
#### Figure 3

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   3A     |   | Barplot of occurrence and co-occurrence of nonsynonymous mutations for the 50 most commonly mutated genes across all tumor types         |
|   3B     |  | Heatmap of co-occurrence and mutual exclusivity of nonsynonymous mutations between genes           |
|   3C     |  | Scatterplot showing correlation of SV and CNV breaks            |
|   3D     |  | Barplot showing chromothripsis frequency for all cancer groups with N >= 3 tumors           |
|   3E     |  | Sina plots of RefSig signature weights across cancer groups         |
#### Figure 4

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   4A     | | ROC curve for TP53 classifier run on stranded RNA-Seq data     |
|   4B     | | Violin/strip plots of TP53 scores across TP53 alteration status           |
|   4C     | |  Violin/strip plots of TP53 expression [log(FPKM)] across TP53 alteration status            |
|   4D     | | Box/strip plots of TP53 scores and telomerase scores across cancer groups          |
|   4E    |  |  Heatmap of RefSig signatures in patients with a hypermutator phenotype       |
|   4F     |  | Forest plot from survival analysis exploring TP53 and telomerase score effects          |
|   4G     |  | Forest plot from survival analysis exploring HGG molecular subtype effects          |
|   4H     |  | Kaplan-Meier cure of HGG tumors by molecular subtype           |

#### Figure 5

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   5A     |  | UMAP of tumors colored by broad histology           |
|   5B     |  | Heatmap of GSVA scores across cancer groups           |
|   5C     |  | Box/strip plots of quanTIseq scores across cancer groups           |
|   5D     |  | Forest plot from survival analysis exploring CD274 expression and immune cell proportion effects           |


### Supplemental Figure Data

#### Figure S2

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   S2A     |             |             |
|   S2B     |             |             |
|   S2C     |             |             |
|   S3D     |             |             |
|   S2E     |             |             |
|   S2F     |             |             |
|   S2G     |             |             |
|   S2H     |             |             |
|   S2I     |             |             |

#### Figure S3

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   S3A     |             |             |
|   S3B     |             |             |
|   S3C     |             |             |
|   S3D     |             |             |
|   S3E     |             |             |


#### Figure S4

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   S4A     |             |             |
|   S4B     |             |             |


#### Figure S5

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   S5A     |             |             |
|   S5B     |             |             |
|   S5C     |             |             |



#### Figure S6

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   S6A     |             |             |
|   S6B     |             |             |
|   S6C     |             |             |
|   S6D     |             |             |
|   S6E     |             |             |
|   S6F     |             |             |


#### Figure S7

| Figure Panel | Filename | Brief Description | Notes on Tabular Data | 
|--------|-------------|--------------------|-------------------------|
|   S7A     |             |             |
|   S7B     |             |             |
|   S7C     |             |             |
|   S7D     |             |             |
|   S7E     |             |             |
|   S7F     |             |             |
|   S7G     |             |             |
|   S7H     |             |             |
|   S7I     |             |             |


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
The `sample_id` impacted are `7316-161` and `7316-85`.
