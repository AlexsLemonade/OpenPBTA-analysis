# Gene Set Enrichment Analysis

Written by Stephanie J. Spielman to supercede previous analyses in `ssGSEA-hallmark`. There are three primary goals:

1. Score hallmark pathways based on expression data using GSVA analysis 
2. Assess consistency between approaches for GSVA score calculation [**Initial results suggest substantial inconsistency**]
3. Analyze scores for highly significant differences among tumor classifications 


Additional documentation is pending analysis module completion.

## Folder Content

+ `01-conduct-gsea-analysis.Rmd` performs the GSVA analysis using RSEM FPMK expression data. Results are saved in two TSV files: `results/gsea_scores_bimodal.tsv` (scores calculated to be bimodal) and `results/gsea_scores_gaussian.tsv`.
+ `02-exploratory-gsea.Rmd` models and begins exploration of results from GSVA analysis. *Script completion and associated results are currently pending.*

