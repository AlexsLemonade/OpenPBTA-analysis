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

`02-HGG-molecular-subtyping-subset-files.R` is a script written to subset the copy number, gene expression, fusion, mutation, histologies' and GISTIC's broad values files to include only high-grade glioma samples. This script produces the relevant subset files that can be found in the `hgg-subset` directory.

`03-HGG-molecular-subtyping-cnv.Rmd` is a notebook written to prepare the copy number data relevant to HGG molecular subtyping. This notebook produces a cnv results table found at `results/HGG_cleaned_cnv.tsv`.

`04-HGG-molecular-subtyping-mutation.Rmd` is a notebook written to prepare the mutation data relevant to HGG molecular subtyping. This notebook produces a mutation results table found at `results/HGG_cleaned_mutation.tsv`.

`05-HGG-molecular-subtyping-fusion.Rmd` is a notebook written to prepare the fusion data relevant to HGG molecular subtyping. This notebook produces a fusion results table found at `results/HGG_cleaned_fusion.tsv`.

`06-HGG-molecular-subtyping-gene-expression.Rmd` is a notebook written to prepare the gene expression data relevant to HGG molecular subtyping. This notebook produces two expression results table _(one for each selection strategy)_ found at `results/HGG_cleaned_expression.polya.tsv` and `HGG_cleaned_expression.stranded.tsv`.

`07-HGG-molecular-subtyping-combine-table.Rmd` is a notebook written to combine the cleaned copy number, mutation, fusion, and gene expression data (prepared in this module's previous notebooks) into one final table of results. This notebook produces one table with the cleaned data found at `results/HGG_cleaned_all_table.tsv` and a table with the molecular subtype information for each HGG sample at `results/HGG_molecular_subtype.tsv`.

`08-1p19q-codeleted-oligodendrogliomas.Rmd` is a notebook written to identify samples in the OpenPBTA dataset that should be classified as 1p/19q co-deleted oligodendrogliomas.

## Folder structure

The structure of this folder is as follows:

```
├── 01-HGG-molecular-subtyping-defining-lesions.Rmd
├── 01-HGG-molecular-subtyping-defining-lesions.nb.html
├── 02-HGG-molecular-subtyping-subset-files.R
├── 03-HGG-molecular-subtyping-cnv.Rmd
├── 03-HGG-molecular-subtyping-cnv.nb.html
├── 04-HGG-molecular-subtyping-mutation.Rmd
├── 04-HGG-molecular-subtyping-mutation.nb.html
├── 05-HGG-molecular-subtyping-fusion.Rmd
├── 05-HGG-molecular-subtyping-fusion.nb.html
├── 06-HGG-molecular-subtyping-gene-expression.Rmd
├── 06-HGG-molecular-subtyping-gene-expression.nb.html
├── 07-HGG-molecular-subtyping-combine-table.Rmd
├── 07-HGG-molecular-subtyping-combine-table.nb.html
├── 08-1p19q-codeleted-oligodendrogliomas.Rmd
├── 08-1p19q-codeleted-oligodendrogliomas.nb.html
├── README.md
├── hgg-subset
│   ├── hgg_focal_cn.tsv.gz
│   ├── hgg_fusion.tsv
│   ├── hgg_gistic_broad_values.tsv
│   ├── hgg_snv_maf.tsv.gz
│   ├── hgg_zscored_expression.polya.RDS
│   └── hgg_zscored_expression.stranded.RDS
├── results
│   ├── HGG_cleaned_all_table.tsv
│   ├── HGG_cleaned_cnv.tsv
│   ├── HGG_cleaned_expression.polya.tsv
│   ├── HGG_cleaned_expression.stranded.tsv
│   ├── HGG_cleaned_fusion.tsv
│   ├── HGG_cleaned_mutation.tsv
│   ├── HGG_defining_lesions.tsv
│   └── HGG_molecular_subtype.tsv
└── run-molecular-subtyping-HGG.sh
```
