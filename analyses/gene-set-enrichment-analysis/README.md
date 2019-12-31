# Gene Set Enrichment Analysis

Written by Stephanie J. Spielman to supercede previous analyses in [`ssgsea-hallmark`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/ssgsea-hallmark). Primary goals include:

1. Score hallmark pathways based on expression data using GSVA analysis, using a strategy that produces Gaussian-distributed scores.
2. Analyze scores for highly significant differences among tumor classifications 


## Folder Content

+ `01-conduct-gsea-analysis.Rmd` performs the GSVA analysis using RSEM FPKM expression data. Results are saved the TSV file `results/gsea_scores.tsv`.
