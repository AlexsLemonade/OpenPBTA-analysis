## OpenPBTA Oncoprint Landscape

**Module authors:** Chante Bethell ([@cbethell](https://github.com/cbethell)) and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

The purpose of this module is to plot the landscape of the consensus CNV, SNV and fusion OpenPBTA datasets.

### Usage

To run the Rscripts in this module from the command line as intended, use:

```
bash run-oncoprint.sh
```

`run-oncoprint.sh` is designed to be run as if it was called from this module directory even when called from outside of this directory.

### Folder content

* `00-map-to-sample_id.R` prepares MAF, focal CN (the output of the `focal-cn-file-preparation` module), and standardized fusion files for use with `01-plot-oncoprint.R`. 
  * The `Tumor_Sample_Barcode` column in the output corresponds to the `sample_id` column in the histologies file
  * We remove ambiguous `sample_id` -- i.e., where there are more than two tumor biospecimens that map to the same sample identifier.
  * Filtering via an [independent specimen file](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples) is optional, but highly recommended.
* `01-plot-oncoprint.R` takes the files from above and optionally a list of top mutated genes (generated in `analyses/interaction-plots/scripts/01-disease-specimen-lists.R`, top genes with recurrent CNVs (generated in `analyses/focal-cn-file-preparation/06-find-recurrent-calls.Rmd`, or a concatenated list of these genes.
This script generates an oncoprint displaying the landscape across PBTA given
the relevant metadata.

### Folder Structure

```
oncoprint-landscape
├── 00-map-to-sample_id.R
├── 01-plot-oncoprint.R
├── README.md
├── data
├── plots
│   ├── all_participants_primary-plus_goi_oncoprint.png
│   ├── all_participants_primary-plus_oncoprint.png
│   ├── all_participants_primary_only_goi_oncoprint.png
│   └── all_participants_primary_only_oncoprint.png
├── run-oncoprint.sh
└── util
    └── oncoplot-functions.R
```
