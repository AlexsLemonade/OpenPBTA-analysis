# Molecular Subtyping non-MB, non-ATRT Embryonal Tumors

**Note: The files in the `subset-files` directory were generated via `02-generate-subset-files.R` using the the files in the [version 13 data release](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/444).**

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Usage](#usage)
- [Folder Content](#folder-content)
  - [Define samples for subsetting](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-embryonal/01-samples-to-subset.nb.html)
  - [Generate subset files](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/molecular-subtyping-embryonal/02-generate-subset-files.R)
  - [Clean C19MC miRNA cluster data](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-embryonal/03-clean-c19mc-data.nb.html)
  - [Construct summary tables](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-embryonal/04-table-prep.nb.html)
- [Folder Structure](#folder-structure)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Usage

When re-running this module, you may want to regenerate the subset files using the most recent data release.
Files will be regenerated using the symlinked files in `data` by default when running from the command line as follows:

```
bash run-embryonal-subtyping.sh
```


## Folder Content

[`01-samples-to-subset.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-embryonal/01-samples-to-subset.nb.html) is a notebook written to identify samples to include in subset files for the purpose of molecularly subtyping non-MB and non-ATRT embryonal tumors.
The samples are identified using the following criteria:

1. An RNA-seq biospecimen sample includes a _TTYH1_ fusion (5' partner).
2. Any sample with "Embryonal tumor" in the `broad_histology` column of the metadata `pbta-histologies.tsv` that is not labeled "Medulloblastoma" or "Atypical Teratoid Rhabdoid Tumor (ATRT)" in `disease_type_old` or `disease_type_new` columns [per this comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/251#issuecomment-568220913).
The output of this notebook is a tsv file containing the biospeciemen IDs identified based on the above criteria (stored in the `results` directory of this module.

[`02-generate-subset-files.R`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/molecular-subtyping-embryonal/02-generate-subset-files.R) is a script written to subset the files required for the subtyping of non-MB and non-ATRT embryonal tumors.
The output of this script includes the subsets of the structural variant, poly-A RNA-seq, and stranded RNA-seq data files (stored in the `subset-files` directory of this module).

[`03-clean-c19mc-data.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-embryonal/03-clean-c19mc-data.nb.html) is a notebook written to clean copy number data related to C19MC amplifications in non-MB, non-ATRT embryonal tumors.
Specifically, the goal of this notebook is to identify embryonal tumors with multilayered rosettes (ETMR), C19MC-altered tumors.
The output of this notebook is a cleaned tsv file containing a binary column indicating whether or not each biospecimen ID is associated with chromosome 19 amplification, which is a suggested characteristic of ETMRs (stored in the `results` directory of this module).

[`04-table-prep.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-embryonal/04-table-prep.nb.html) is a notebook written to construct tables that summarize the data relevant to the molecular subtyping of non-MB and non-ATRT embryonal tumors per the information provided on the [reference GitHub issue](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/251).
The output of this notebook includes two tsv files, both found in the `results` directory of this module.
The first output file, `embryonal_tumor_subtyping_relevant_data.tsv`, contains a summary table of the data subsetted in `02-generate-subset-files.R` for the purpose of molecular subtyping embryonal tumors.  
The second output file, `embryonal_tumor_molecular_subtypes.tsv`, contains the molecular subtype information of the identified biospecimen IDs based on the summarized relevant data as described in the [original comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/251#issue-520154478) and in [this comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/251#issuecomment-571807158) both found on the reference GitHub issue.
The information in this file is represented in table with the following columns:


| `Kids_First_Participant_ID` | `sample_id` | `Kids_First_Biospecimen_ID_DNA` | `Kids_First_Biospecimen_ID_RNA` | `molecular_subtype` |
|-----------------------------|-------------|---------------------------------|---------------------------------|---------------------|



## Folder Structure

```
├── 01-samples-to-subset.Rmd
├── 01-samples-to-subset.nb.html
├── 02-generate-subset-files.R
├── 03-clean-c19mc-data.Rmd
├── 03-clean-c19mc-data.nb.html
├── 04-table-prep.Rmd
├── 04-table-prep.nb.html
├── README.md
├── results
│   ├── biospecimen_ids_embryonal_subtyping.tsv
│   ├── cleaned_chr19_cn.tsv
│   ├── embryonal_tumor_molecular_subtypes.tsv
│   └── embryonal_tumor_subtyping_relevant_data.tsv
├── run-embryonal-subtyping.sh
└── subset-files
    ├── embryonal_manta_sv.tsv
    ├── embryonal_zscored_exp.polya.rds
    └── embryonal_zscored_exp.stranded.rds
```
