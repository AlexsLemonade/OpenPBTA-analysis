# Mutational Signatures

This analysis evaluates mutational signatures of the [consensus SNV callers file](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers#consensus-mutation-call).

The overall analysis considers three different approaches to assessing mutational signatures in PBTA data:

+ Evaluate signatures from [COSMIC signatures](https://cancer.sanger.ac.uk/cosmic)
and [Alexandrov et al, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23945592) for all samples using [`deconstructSigs`](https://github.com/raerose01/deconstructSigs).

+ Conduct `de novo` signature extraction from WGS samples using [`sigfit`](https://github.com/kgori/sigfit), including goodness-of-fit analyses to determine optimal parameters.
  + This analysis requires up to 16 GB RAM.

+ Evaluate known (adult) CNS signatures, as identified in [Degasperi et al, 2020](https://doi.org/10.1038/s43018-020-0027-5) in WGS samples using [`sigfit`](https://github.com/kgori/sigfit).


## Usage

To reproduce this analysis, run the following command
```
bash run-mutational-signatures.sh
```

## Contents

+ The Rmd `01-known_signatures.Rmd` runs the mutational signatures analysis using existing COSMIC and Alexander et al, 2013 signatures. 
  + Calculates weights from [`deconstructSigs::whichSignatures`](https://www.rdocumentation.org/packages/deconstructSigs/versions/1.8.0/topics/whichSignatures) which are are multiplied by each sample's sum of mutations as provided by [`deconstructSigs::mut.to.sigs.input`](https://www.rdocumentation.org/packages/deconstructSigs/versions/1.8.0/topics/mut.to.sigs.input). These numbers are saved to the `_signatures_results.tsv` files in the `results` folder.

+ The script `02-split_experimental_strategy.R` splits up the consensus MAF files by experimental strategy, and writes WGS_only data to `scratch/` for consumption by the rest of the mutational signature analysis

+ The script `03-de_novo_range_of_nsignatures.sh` is run twice. First, it is used perform benchmarking and ascertain appropriate `k` (number of signatures to extract) and robustness to random seed for *de novo* signature extraction
  + The benchmarking procedure exports several PNG images containing goodness-of-fit cosine plots to `plots/denovo/gof/` which are visually used to assess optimal `k` values. These files are named according to random seed and the model used for extraction (`multinomial` which corresponds to Alexander et al 2013 or `poisson` which corresponds to the approach from [Emu](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-4-r39)).
  
+ The script `03-de_novo_range_of_nsignatures.sh`  is again run to perform the *de novo* signature extraction with the determined range of `k` values.
  + The extraction procedure exports both PNG images containing (refined) goodness-of-fit cosine plots to `plots/denovo/extraction/`and full fitted results containing fitted signatures (`results/de_novo_signatures.RDS`) and associated exposures (`results/de_novo_exposures.RDS`).

+ The Rmd `04-analyze_de_novo.Rmd` analyzes the results from *de novo* signature extraction. This file knits to an HTML (`04-analyze_de_novo.HTML`) which contains a **table** and description of extraction signatures.

+ The script `05-fit_cns_signatures.R` determines the relative contributions of 8 known (adult) CNS signatures in the PBTA WGS data.
  + This script saves the fitted CNS signature exposures in `results/fitted_cns_signature_exposures.RDS`.

+ The Rmd `06-analyze_cns_fit.Rmd` analyzes and visualizes the results from CNS signature fitting. 
  + The knitted HTML contains several figures which are also separately exported to `plots/cns/`:
    + `mean_median_exposure_barplot.png` displays the mean and median exposures for CNS signatures for each of the main histology groups we consider in PBTA.
    + `top_10_samples_barplot.png` displays the relative signature exposures for the top-10 most mutated samples in the PBTA.


## Overall file structure

```
OpenPBTA-analysis
├── analyses
│   └── mutational-signatures
│       ├── mutational_signatures.Rmd
│       ├── util
│       │    └── mut_sig_functions.R
│       ├── results
│       │   ├── cosmic_signatures_results.tsv
│       │   └── nature_signatures_results.tsv
│       │   └── de_novo_exposures.RDS
│       │   └── de_novo_signatures.RDS
│       │   └── fitted_cns_signature_exposures.RDS
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
│           |   └── ...
│           └── denovo
│           │   ├── gof
│           │   │   ├── seed_<RANDOM SEED>_model_<EXTRACTION MODEL>.png
│           │   │   └── ...
│           │   ├── extraction
│           │   │   ├── seed_<RANDOM SEED>_model_<EXTRACTION MODEL>.png
│           │   │   └── ...
│           └── cns
│           │   ├── mean_median_exposure_barplot.png
│           │   └── top_10_samples_barplot.png

├── data
```

## Summary of custom functions

|Function Name|Summary|
|-------------|-----------|
|`sample_mut_sig_plot`|Saves traditional mutational signature plots for each sample|
|`calc_mut_per_sig`|Given `deconstructSigs::whichSignature` output, formats the sample data into a data.frame and calculates the mutations per Mb for each sample and each signature|
|`bubble_matrix_plot`|Groups together data by histology and makes the bubble matrix plot|
|`grouped_sig_barplot`|Creates signature grouped barplots for histology group provided and only uses primary tumors' data|
