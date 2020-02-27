# GISTIC Comparision of Entire Cohort Results vs Individual Histology Results

**Module authors:** Chante Bethell ([@cbethell](https://github.com/cbethell)) and Jaclyn Taroni ([@jaclyn-taroni](https://github.com/jaclyn-taroni))

**Note: The files used to run the notebook in this directory were generated via [GISTIC](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/doc/data-formats.md#gistic-output-file-formats) and stored in the `analyses/run-gistic/results` directory of this repository.
When re-running this module, you may want to regenerate the GISTIC output files using the most recent data release.**

## Usage

To run the R notebook in this module from the command line (assuming that you are in the top directory of this repository), use:

```
Rscript -e "rmarkdown::render('analyses/gistic-cohort-vs-histology-comparison/gistic-cohort-vs-histology-comparison.Rmd', clean = TRUE)"
```

## Folder content

[`gistic-cohort-vs-histology-comparison.nb.html`](https://alexslemonade.github.io/OpenPBTA-analysis/analyses/gistic-cohort-vs-histology-comparison/gistic-cohort-vs-histology-comparison.nb.html) is a notebook written to identify, if any, disagreement between GISTIC results for the entire cohort versus the three individual histologies that we have GISTIC results for (LGAT, HGAT, and medulloblastoma).

## Folder Structure

```
├── README.md
├── gistic-cohort-vs-histology-comparison.Rmd
├── gistic-cohort-vs-histology-comparison.nb.html
└── plots
    ├── hgat_gistic_scores_multipanel_plot.png
    ├── lgat_gistic_scores_multipanel_plot.png
    └── medulloblastoma_gistic_scores_multipanel_plot.png
```
