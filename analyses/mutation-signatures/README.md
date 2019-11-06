# Mutation Signatures

This analysis evaluates mutation signatures of the [consensus SNV callers file](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers#consensus-mutation-call).

Here the signatures from [COSMIC signatures](https://cancer.sanger.ac.uk/cosmic) and [Alexandrov et al, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23945592) are evaluated for all samples using [deconstructSigs](https://github.com/raerose01/deconstructSigs).

## Summary of Results:

Coming soon.

**Table of Contents**
* [How to run this analysis](#how-to-run-the-mutation-signatures-analysis)
* [Summary of the calculations](#summary-of-the-calculations)
* [Overall file structure](#overall-file-structure)
* [Summary of functions](#summary-of-custom-functions)

## How to run the mutation signatures analysis

To run the mutation signature evaluations, use this in command line:
```
Rscript -e "rmarkdown::render('analyses/mutation-signatures/mutation_signatures.Rmd',
                              clean = TRUE)"
```

This notebook, in addition to its `nb.html` output, will return plots return results for each caller in the `plots` and `results` folder.

## Summary of the calculations

### Number of mutations per signature

The weights from [`deconstructSigs::whichSignatures`](https://www.rdocumentation.org/packages/deconstructSigs/versions/1.8.0/topics/whichSignatures) are multiplied by each sample's sum of mutations as provided by [`deconstructSigs::mut.to.sigs.input`](https://www.rdocumentation.org/packages/deconstructSigs/versions/1.8.0/topics/mut.to.sigs.input).
These numbers are saved to the `_signatures_results.tsv` files in the `results` folder.
For more information, see the [`calc_mut_per_sig`]() code.

### Proportion of tumors with a signature

A tumor is considered to have a particular signature if its weight is non-zero.
This is divided by the number of tumor samples in that particular histology group.
For more information, see the [`bubble_matrix_plot`]() code.

## Overall file structure
```
OpenPBTA-analysis
├── analyses
│   └── mutation-signatures
│       ├── mutation_signatures.Rmd
│       ├── util
│       │    └── mut_sig_functions.R
│       ├── results
│       │   ├── cosmic_signatures_results.tsv
│       │   └── nature_signatures_results.tsv
│       └── plots
│           ├── cosmic
│           │   ├── individual_mut_sigs
│           │   │   ├── <TUMOR_SAMPLE_BARCODE>_cosmic_mutation_sig.png
│           │   │   └── ...
│           │   ├── signature_grouped_barplots
│           │   │   ├── barplot_<HISTOLOGY_GROUP>_cosmic_mutation_sig.png
│           │   │   └── ...
│           │   └── bubble_matrix_cosmic_mutation_sig.png
│           └── nature
│               └── ...
├── data
```

## Summary of custom functions

|Function Name|Summary|
|-------------|-----------|
|`sample_mut_sig_plot`|Saves traditional mutation signature plots for each sample|
|`calc_mut_per_sig`|Given `deconstructSigs::whichSignature` output, formats the sample_data into a data.frame and calculates the mutations per Mb for each sample and each signature|
|`bubble_matrix_plot`|Groups together data by histology and makes the bubble matrix plot|
|`grouped_sig_barplot`|Creates signature grouped barplots for histology group provided and only uses primary tumors' data|
