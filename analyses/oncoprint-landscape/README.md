## OpenPBTA Oncoprint Landscape

**Module authors:** Chante Bethell ([@cbethell](https://github.com/cbethell)) and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

The purpose of this module is to plot the landscape of the consensus CNV, SNV and fusion OpenPBTA datasets.

### Usage

To run the Rscripts in this module from the command line as intended, use:

```
bash run-oncoprint.sh
```

`run-oncoprint.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

### Individual scripts

* `00-prepare-goi-lists.R` prepares genes of interest (GOI) lists from `data/oncoprint-goi-lists-OpenPBTA.tsv`, as well as a table that maps between cancer groups and the relevant GOI list. 
All output is in `data/` (within this module).
* `01-map-to-sample_id.R` prepares MAF, focal CN (from the "most focal" output of the `focal-cn-file-preparation` module), and standardized fusion files for use with `01-plot-oncoprint.R`.
  * The `Tumor_Sample_Barcode` column in the output corresponds to the `sample_id` column in the histologies file
  * We remove ambiguous `sample_id` -- i.e., where there are more than two tumor biospecimens that map to the same sample identifier.
  * Filtering via an [independent specimen file](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples) is optional, but highly recommended.
* `02-plot-oncoprint.R` takes the files from above and optionally a file or set of files (to be concatenated) that will restrict the set of genes that are being plotted in an OncoPrint.
  * This script also returns tables summary of the total of number of samples that have alterations in gene of interest in `tables/`
* `03-oncoprint-n-count-table.R` counts the number of samples that enter oncoprint plotting (with or without alterations in genes of interest) for each broad histology group.
The output is in `tables/`.
* `04-alteration-counts-by-cancer-group.R` summarizes the alterations input into the oncoprint plotting, when genes of interest (GOI) lists are used, for each cancer group (with non-synonymous mutations after subsetting) separately and outputs this to the `tables/cancer_group_counts/` directory.

_No longer used in any analysis in **this** module_: Gene lists in driver-lists were obtained from https://github.com/marislab/create-pptc-pdx-oncoprints/tree/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data and are still used in `focal-cn-file-preparation`.