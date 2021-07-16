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

* `01-map-to-sample_id.R` prepares MAF, focal CN (from the "most focal" output of the `focal-cn-file-preparation` module), and standardized fusion files for use with `01-plot-oncoprint.R`. 
  * The `Tumor_Sample_Barcode` column in the output corresponds to the `sample_id` column in the histologies file
  * We remove ambiguous `sample_id` -- i.e., where there are more than two tumor biospecimens that map to the same sample identifier.
  * Filtering via an [independent specimen file](https://alexslemonade.github.io/OpenPBTA-manuscript/#selection-of-independent-samples) is optional, but highly recommended.
* `02-plot-oncoprint.R` takes the files from above and optionally a file or set of files (to be concatenated) that will restrict the set of genes that are being plotted in an OncoPrint.
	* Running this via `run-oncoprint.sh` will restrict plotting to a list of top mutated genes (generated in `analyses/interaction-plots/scripts/01-disease-specimen-lists.R`) and top genes with recurrent CNVs (generated in `analyses/focal-cn-file-preparation/06-find-recurrent-calls.Rmd`)


### Folder Structure

```
├── 01-map-to-sample_id.R
├── 02-plot-oncoprint.R
├── README.md
├── driver-lists
│   ├── brain-goi-list-long.txt
│   └── brain-goi-list-short.txt
├── plots
│   ├── primary-plus_goi_oncoprint.png
│   ├── primary-plus_oncoprint.png
│   ├── primary_only_goi_oncoprint.png
│   └── primary_only_oncoprint.png
├── run-oncoprint.sh
└── util
    └── oncoplot-functions.R
```

Gene lists in `driver-lists` were obtained from https://github.com/marislab/create-pptc-pdx-oncoprints/tree/2c9ed2a2331bcef3003d6aa25a130485d76a3535/data.
