# Tumor Mutation Burden Compare to TCGA

This analysis comparees the Pediatric Brain Tumor samples of this dataset to the adult brain tumor samples from TCGA. 

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

### PBTA Tumor Mutation Burden

### TCGA Tumor Mutation Burden
