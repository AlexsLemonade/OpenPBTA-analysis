# Molecular Subtyping High-grade Gliomas

**Module authors:** Chante Bethell ([@cbethell](https://github.com/cbethell)), Stephanie J. Spielman([@sjspielman](https://github.com/sjspielman)), and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

**Note: The files in the `hgg-subset` directory were generated via `02-HGG-molecular-subtyping-subset-files.R` using the the files in the [version 13 data release](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/444).
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

[`01-HGG-molecular-subtyping-defining-lesions.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-HGG/01-HGG-molecular-subtyping-defining-lesions.nb.html) is a notebook written to look at the high-grade glioma defining lesions for all samples in the PBTA dataset and isolates samples that have been classified as high-grade gliomas and that may be reclassified as high-grade gliomas based on having mutations at at least one of the defining lesions.
This notebook produces a results table found at `results/HGG_defining_lesions.tsv`.

`02-HGG-molecular-subtyping-subset-files.R` is a script written to subset the copy number, gene expression, fusion, mutation, SNV and GISTIC's broad values files to include only high-grade glioma samples.
This script produces the relevant subset files that can be found in the `hgg-subset` directory.

[`03-HGG-molecular-subtyping-cnv.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-HGG/03-HGG-molecular-subtyping-cnv.nb.html) is a notebook written to prepare the copy number data relevant to HGG molecular subtyping.
The CNVkit focal copy number file generated in the [`focal-cn-file-preparation`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/focal-cn-file-preparation/) module is used as this CNVkit was also used to produce the GISTIC `broad_values_by_arm.txt` file that is also implemented in this module.
The decision was also made to use consistent coding of the values in the output table as follows: `"loss"`, `"neutral"`, and `"gain"`.
This notebook produces a cnv results table with cleaned CNVkit and GISTIC data found at `results/HGG_cleaned_cnv.tsv`.

[`04-HGG-molecular-subtyping-mutation.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-HGG/04-HGG-molecular-subtyping-mutation.nb.html) is a notebook written to prepare the consensus mutation data relevant to HGG molecular subtyping.
Per [issue #249](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/249) in the OpenPBTA-analysis repository, we filtered the subsetted SNV data to the genes of interest based on the following criteria: 1) for the `H3 K28 mutant` subtype, cooccurring lesions include _ACVR1, TP53, ATRX_ mutations, mutually exclusive lesions include _FGFR1_ and _IDH1_ mutations and _BRAF V600E_ (specific for low-grade gliomas);
2) for the `H3 G35 mutant` subtype, cooccurring lesions include _ATRX/DAXX, TP53, SETD2_ mutations and mutually exclusive lesions include _IDH1_ mutations;
3) for the `IDH mutant` subtype,  _IDH1_ R132H mutations are harbored and cooccurring lesions include _TP53_ mutations;
4) for the `H3.3 and IDH wildtype` subtype, _H3F3A_ or _IDH1_ mutations are absent and _TP53_ and _TERT_ mutations are harbored;
5) for the `1p/19q co-deleted oligodendrogliomas` subtype, nearly all have a cooccurring _IDH1_ mutation (R132 or R172) and other cooccurring lesions include _TERT_ promoter, _CIC_, and _FUBP1_ mutations.
This notebook wrangles the SNV data and produces a consensus mutation results table, representing the above information in a way that can be used to molecularly subtype the HGG samples.
This table is found at `results/HGG_cleaned_mutation.tsv`.

[`05-HGG-molecular-subtyping-fusion.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-HGG/05-HGG-molecular-subtyping-fusion.nb.html) is a notebook written to prepare the putative oncogenic fusion data relevant to HGG molecular subtyping.
Per [issue #249](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/249), we filtered the data to the two fusions of interest: _FGFR1_ fusions, which should be mutually exclusive of H3 K28 mutants, and _NTRK_ fusions, which are co-occurring with H3 G35 mutants.
There is no mention of specific fusion partners or orientations, so we look at _any_ instances of fusions that include _FGFR1_ or _NTRK_.
_Note:_ _NTRK_ refers to a [family of receptor kinases](https://www.biooncology.com/pathways/cancer-tumor-targets/ntrk/ntrk-oncogenesis.html), so we include the full fusion name to account for various individual _NTRK_ gene symbols.
This notebook produces a fusion results table found at `results/HGG_cleaned_fusion.tsv`.

[`06-HGG-molecular-subtyping-gene-expression.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-HGG/06-HGG-molecular-subtyping-gene-expression.nb.html) is a notebook written to prepare the gene expression data relevant to HGG molecular subtyping.
Per [issue #249](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/249), we filtered the gene expression data to three genes of interest: _OLIG2_ and _FOXG1_ should be highly expressed in IDH mutants, and _TP73-AS1_ methylation and downregulation cooccurs with _TP53_ mutations.
This notebook produces two expression results table _(one for each selection strategy)_ found at `results/HGG_cleaned_expression.polya.tsv` and `HGG_cleaned_expression.stranded.tsv`.

[`07-HGG-molecular-subtyping-combine-table.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-HGG/07-HGG-molecular-subtyping-combine-table.nb.html) is a notebook written to combine the cleaned copy number, mutation, fusion, and gene expression data (prepared in this module's previous notebooks) into one final table of results.
This notebook produces one table with the cleaned data found at `results/HGG_cleaned_all_table.tsv`.
A table with the molecular subtype information for each HGG sample at `results/HGG_molecular_subtype.tsv` is also produced, where the subtype values in the `molecular_subtype` column are determined as follows:
1) If there was an _H3F3A_ K28M or _HIST1H3B_ K28M mutation and no _BRAF_ V600E mutation -> `HGG_H3_K28`;
2) If there was an _H3F3A_ K28M or _HIST1H3B_ K28M mutation and a _BRAF_ V600E mutation -> `BRAF_V600E_mutation`;
3) If there was an _H3F3A_ G35V or G35R mutation -> `HGG_H3_G35`;
4) If there was an _IDH1_ R132 mutation -> `HGG_IDH`.

[`08-1p19q-codeleted-oligodendrogliomas.Rmd`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/molecular-subtyping-HGG/08-1p19q-codeleted-oligodendrogliomas.nb.html) is a notebook written to identify samples in the OpenPBTA dataset that should be classified as 1p/19q co-deleted oligodendrogliomas.
The GISTIC `broad_values_by_arm.txt` file is used to identify samples with `1p` and `19q` loss, then the consensus mutation file is filtered to the identified samples in order to check for _IDH1_ mutations.
_Note:_ Per [this comment](https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/435#issuecomment-576898275), very few samples in the OpenPBTA dataset, if any, are expected to fit into the `1p/19q co-deleted oligodendrogliomas` subtype.

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
