# ssGSEA Hallmark Updated

Written by Stephanie J. Spielman to supercede previous analyses in `ssGSEA-hallmark`. There are two primary goals:

1. Score hallmark pathways based on expression data using ssGSEA analysis (complete pending code review)
2. Analyze scores for highly significant differences among tumor classifications (20% complete)


Additional documentation is pending analysis module completion.

## Folder Content

+ `01-conduct-ssGSEA-analysis.Rmd` performs the ssGSEA analysis using RSEM FPMK expression data. Results are saved in a TSV `results/ssGSEA_scores.tsv`
+ `02-model-ssGSEA-results.Rmd` models and visualizes results from the ssGSEA analysis. *Script completion and associated results are currently pending.*
