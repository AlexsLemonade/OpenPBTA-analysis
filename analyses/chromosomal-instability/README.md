# Chromosomal Instability Analysis

This analysis evaluates chromosomal instability by formatting, overlapping and
then mapping SV and CNV data.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Usage](#usage)
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

## Output

This notebook returns chromosomal break plots for each sample in the `plots` directory.

You can see the output notebook [here](https://cansavvy.github.io/openpbta-notebook-concept/chromosomal-instability/chromosomal-instability.nb.html).

## Summary of Custom Functions

- `make_granges`
- `break_density`
- `map_density_plot`
