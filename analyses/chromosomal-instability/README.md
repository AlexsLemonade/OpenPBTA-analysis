# Chromosomal Instability Analysis

This analysis evaluates chromosomal instability by calculating chromosomal
break point density calculations and circular plot visualization using CNV and SV data.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Usage](#usage)
- [Methods](#methods)
- [Output](#output)
- [Summary of Custom Functions](#summary-of-custom-functions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Usage

This analysis can be run via the command line from the top directory of the
repository as follows:

```
bash run-breakpoint-analysis.sh
```

## Methods

CNV and SV data are used to calculate chromosomal instability.
Then both CNV and SV range datasets are transformed into single breakpoint data.
Breakpoint density is calculated by creating bins using `GenomicRanges::tileGenome` using a one Mb window size.
The `01-plot-chromosomal-instability.Rmd` returns chromosomal break plots for each sample and `short_histology` group in the `plots` directory.
`01b-visualization-cnv-sv.Rmd` notebook shows examples of how to plot CNV and SV data to visualize locations of chromosomal instability.

## Output

Three output TSVs (one for each the CNV and SV data, and one for the combined data) with breakpoint density per Mb of the effectively surveyed genome are saved to `breakpoint-data` directory.
Note that this script is set up to handle WXS and WGS separately, however currently our PBTA dataset only has CNV and SV data for WGS samples.
The individual sample plots and grouped by `short_histology` plots are in the `plots/sample` and `plots/tumor-type` directories, respectively.

## Summary of Custom Functions

For breakpoint analysis:
- `make_granges` : Given a data.frame with chr break coordinates, make a `GenomicRanges` object.
- `break_density`: Given data.frame(s) with chr break coordinates, calculate the density of the breaks.
- `map_breaks_plot`: Given a `GenomicRanges` object, use map the chromosomal coordinates to a `ggplot2`
- `multipanel_break_plot`: Given a list of `GenomicRanges` objects, plot them in a combined `cowplot`.

For circular plotting:
`circos_map_plot`: Given a data.frame with chromosomal coordinates, and a corresponding data value to plot, make a circos plot or add a circos track to an existing plot.
`circos_map_transloc`: Given a data.frame with two sets of coordinates, map the links between those coordinates on a new or existing circos plot.
`prep_bed`: Internally used by the circos functions to make sure data is prepped for `circlize`.
