# Chromosomal Instability Analysis

This analysis evaluates chromosomal instability by calculating chromosomal
break point densities using CNV and SV data.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Usage](#usage)
- [Methods](#methods)
- [Output](#output)
- [Summary of Custom Functions](#summary-of-custom-functions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Usage

This notebook can be run via the command line from the top directory of the
repository as follows:

```
Rscript -e "rmarkdown::render('analyses/chromosomal-instability/chromosomal-instability.Rmd',
                              clean = TRUE)"
```

## Methods

CNV and SV data are used to calculate chromosomal instability.
SV data here is first formatted by running `sv-analysis/01-process-sv-file.R`
script and using the without Y and M files.
Then both CNV and SV range datasets are transformed into single breakpoint data.
Breakpoint density is calculated by creating bins using `GenomicRanges::tileGenome` using a one Mb window size.
This notebook returns chromosomal break plots for each sample and `short_histology` group in the `plots` directory.

## Output

The individual sample plots and grouped by `short_histology` plots are in the `plots/sample` and `plots/tumor-type` directories, respectively.

## Summary of Custom Functions

- `make_granges` : Given a data.frame with chr break coordinates, make a `GenomicRanges` object.
- `break_density`: Given data.frame(s) with chr break coordinates, calculate the density of the breaks.
- `map_density_plot`: Given a `GenomicRanges` object, use map the chromosomal coordinates to a `ggplot2`
- `chr_break_plot`: Given a list of `GenomicRanges` objects, plot them in a combined `cowplot`.
