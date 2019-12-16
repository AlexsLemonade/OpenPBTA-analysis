# Tumor Mutation Burden Compare to TCGA

This analysis compares the Pediatric Brain Tumor samples of this dataset to the adult brain tumor samples from [The Cancer Genome Atlas](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Usage](#usage)
- [Results](#results)
- [Summary of Methods](#summary-of-methods)
  - [PBTA Tumor Mutation Burden](#pbta-tumor-mutation-burden)
  - [TCGA Tumor Mutation Burden](#tcga-tumor-mutation-burden)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Usage

To run this from the command line, use:
```
Rscript -e "rmarkdown::render('analyses/tmb-compare-tcga/compare-tmb.Rmd',
                              clean = TRUE)"
```

## Results

![](plots/tmb_tcga_and_pbta_plot.png)

## Summary of Methods

In brief, TMB for both brain tumor datasets are calculated using mutations from exonic regions of the genome only.  

### PBTA Tumor Mutation Burden

The TMB calculations for the pediatric brain tumor set were carried out in [snv-callers analysis](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) in this repository.
In brief, tumor mutation burden is calculated using consensus SNV calls found in exonic regions.
Consensus mutations included any calls that were made by all three of: [Mutect2](https://software.broadinstitute.org/cancer/cga/mutect), [Strelka2](https://github.com/Illumina/strelka), and [Lancet](https://github.com/nygenome/lancet).
The total number exonic consensus SNVs were used for the numerator and the effective size of the genome surveyed is used as the denominator.

```
TMB = (total # exonic consensus snvs) / (size of effectively surveyed genome)
```

For more details, see [snv-callers README](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers).

### TCGA Tumor Mutation Burden

For calculating TCGA tumor mutation burden, [MC3 mutation calls](https://gdc.cancer.gov/about-data/publications/mc3-2017) were used with TCGA Brain tumor related projects only:

- [LGG (Lower-grade Glioma)](https://www.nejm.org/doi/full/10.1056/NEJMoa1402121)
- [GBM (Glioblastoma Multiforme)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga/studied-cancers/glioblastoma)
- [PCPG (Pheochromocytoma and Paraganglioma)]()

The size of the [exome BED regions file](https://api.gdc.cancer.gov/data/7f0d3ab9-8bef-4e3b-928a-6090caae885b) included with the MC3 project was used for the denominator.
All mutations for all three projects fell within these exome bed regions.
