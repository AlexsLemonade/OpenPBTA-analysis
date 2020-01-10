# Molecular Subtyping High-grade Gliomas

**Note: The files in the `hgg-subset` directory were generated via `02-HGG-molecular-subtyping-subset-files.R` using the the files in the [version 12 data release](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/347).
When re-running this module, you may want to regenerate the HGG subset files using the most recent data release.**

## Usage

To run all of the Rscripts in this module from the command line sequentially, use:

```
bash run-molecular-subtyping-HGG.sh
```

When run in this manner, `02-HGG-molecular-subtyping-subset-files.R` will generate subset files using whichever files are symlinked in `data` on your local machine.

`run-molecular-subtyping-HGG.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

## Folder content

This folder contains scripts tasked to molecularly subtype High-grade Glioma samples in the PBTA dataset.

`01-HGG-molecular-subtyping-defining-lesions.Rmd` is a notebook written to look at the high-grade glioma defining lesions for all samples in the PBTA dataset and isolates samples that have been classified as high-grade gliomas and that may be reclassified as high-grade gliomas based on having mutations at at least one of the defining lesions. This notebook produces a results table found at `results/HGG_defining_lesions.tsv`.

`02-HGG-molecular-subtyping-subset-files.R` is a script written to subsets the focal copy number, RNA expression, fusion, histologies' and GISTIC's broad values files to include only high-grade glioma samples. This script produces the relevant subset files that can be found in the `hgg-subset` directory.

`03-HGG-molecular-subtyping-data-prep.Rmd` is a notebook written to prepare the copy number, RNA expression, structural variant, fusion, and GISTIC broad values by chromosomal arm data that are relevant for molecular subtyping. This notebook produces a final results table found at `results/HGG_molecular_subtypes.tsv`. 

## Folder structure

The structure of this folder is as follows:

```
├── 01-HGG-molecular-subtyping-defining-lesions.Rmd
├── 01-HGG-molecular-subtyping-defining-lesions.nb.html
├── 02-HGG-molecular-subtyping-subset-files.R
├── 03-HGG-molecular-subtyping-data-prep.Rmd
├── 03-HGG-molecular-subtyping-data-prep.nb.html
├── README.md
├── hgg-subset
│   ├── hgg_focal_cn.tsv.gz
│   ├── hgg_fusion.tsv
│   ├── hgg_gistic_broad_values.tsv
│   ├── hgg_histologies.tsv
│   └── hgg_zscored_expression.RDS
├── results
│   ├── HGG_defining_lesions.tsv
│   └── HGG_molecular_subtypes.tsv
└── run-molecular-subtyping-HGG.sh
```
