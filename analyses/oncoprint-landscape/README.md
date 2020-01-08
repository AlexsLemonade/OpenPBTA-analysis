# OpenPBTA Oncoprint

## Work in progress

This module is a work in progress.
The code was written to be updated as consensus, filtered, and/or prioritized data becomes available.
**The plots included as PNGs should be regarded as proof of concept, rather than for interpretation.**

In particular, the following are likely to change:

* The plot will use consensus copy number calls ([#128](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/128)).

## Usage

To run the Rscript in this module from the command line as intended, use:

```
bash run-oncoprint.sh
```

`run-oncoprint.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

## Folder content

* `00-map-to-sample_id.R` prepares MAF, focal CN (the output of the `focal-cn-file-preparation` module), and standardized fusion files for use with `01-plot-oncoprint.R`. 
  * The `Tumor_Sample_Barcode` column in the output corresponds to the `sample_id` column in the histologies file
  * We remove ambiguous `sample_id` -- i.e., where there are more than two tumor biospecimens that map to the same sample identifier.
  * Filtering via an [independent specimen file](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples) is optional, but highly recommended.
* `01-plot-oncoprint.R` takes the files from above and optionally a gene list and creates an oncoprint.

