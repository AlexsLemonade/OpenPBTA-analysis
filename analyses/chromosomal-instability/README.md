# Chromosomal Instability Analysis

This analysis evaluates chromosomal instability by calculating chromosomal
break point density calculations and circular plot visualization using CNV and SV data.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Usage](#usage)
- [Methods and Output](#methods-and-output)
- [Summary of Custom Functions](#summary-of-custom-functions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Usage

This analysis can be run via the command line from the top directory of the
repository as follows:

```
bash analyses/chromosome-instability/run_breakpoint_analysis.sh
```

## Methods and Output

- `00-setup-breakpoint-data.R` CNV and SV data are transformed into single breakpoint data as well as their intersection.
Manta SV calls are filtered for PASS variants.
These three datasets (intersection, CNV, and SV) are used to calculate genome wide breakpoint densities for each sample which are saved to three `_breaks_densities.tsv` files. Note that this script is set up to handle WXS and WGS separately, however **currently our PBTA dataset only has CNV and SV data for WGS samples**.
Additionally, a  `breakpoint-data/breaks_lists.RDS` file is saved which contains three data.frames:
1) `intersection_of_breaks` contains the intersection break counts for both SV and CNV break data.  
2) `cnv_breaks` contains the number of break counts for CNV.   
3) `sv_breaks` contains the number of break counts for SV.  

- `01-localization-of-breakpoints.Rmd` uses the data in `breaks_lists.RDS` to co-localize and map breakpoints by bins across the genome.
These binned breakpoint counts are calculated by sample as well as by histology group.  
Bins are created using `GenomicRanges::tileGenome` using a one Mb window size.
Genome bins above a percentage (default is 75%) of their total size being covered in uncallable regions are called as NA for all output statistics.  
The output of this notebook is three `_binned_breakpoint_counts.tsv"` for each dataset, and an RDS file with the histology binned data: `histology_breakpoint_densities.RDS`.

- `02a-plot-chr-instability-heatmaps.Rmd` uses `_binned_breakpoint_counts.tsv"` datasets to create three heatmaps for the intersection, CNV, and SV data respectively. `NA` regions are gray.

- `02b-plot-chr-instability-by-histology.Rmd` uses the `_breaks_densities.tsv` files and `histology_breakpoint_binned_counts.RDS` to plot breakpoint densities by `short_histology` group. 

## Summary of Custom Functions

For breakpoint analysis:
- `make_granges` : Given a data.frame with chr break coordinates, make a `GenomicRanges` object.
- `break_density`: Given data.frame(s) with chr break coordinates, calculate the density of the breaks.
- `map_breaks_plot`: Given a `GenomicRanges` object, use map the chromosomal coordinates to a `ggplot2`
- `multipanel_break_plot`: Given a list of `GenomicRanges` objects, plot them in a combined `cowplot`.
- `breaks_cdf_plot`: Given a genome wide breaks density file path, plot the CDF distribution for it by histology.
