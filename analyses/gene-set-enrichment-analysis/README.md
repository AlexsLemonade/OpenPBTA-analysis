# Gene Set Enrichment Analysis

Written by Stephanie J. Spielman to supercede previous analyses in [`ssgsea-hallmark`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/ssgsea-hallmark). Primary goals include:

1. Score hallmark pathways based on expression data using GSVA analysis, using a strategy that produces Gaussian-distributed scores.
2. Analyze scores for highly significant differences among tumor classifications 
3. [Pending] generate heatmaps and templates for boxplot figures 


## Folder Content

+ `01-conduct-gsea-analysis.R` performs the GSVA analysis using RSEM FPKM expression data for both stranded and polyA data. Results are saved in `results/` TSV files.

+ `02-model-gsea.Rmd` performs ANOVA and Tukey tests on GSVA scores to evaluate, for each hallmark pathway, differences in GSVA across groups (e.g. short histology or disease type).

+ FORTHCOMING: 	`03-visualize-gsea.Rmd` will create heatmaps and boxplots to visualize scores.

+ `results/gsva_scores_stranded.tsv` represents GSVA scores calculated from `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` (data release v12)
	+ File created with: `Rscript --vanilla 01-conduct-gsea-analysis.R --input pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output gsva_scores_stranded.tsv`
+ `results/gsva_scores_polya.tsv` represents GSVA scores calculated from `pbta-gene-expression-rsem-fpkm-collapsed.polya.rds` (data release v12)
	+ File created with: `Rscript --vanilla 01-conduct-gsea-analysis.R --input pbta-gene-expression-rsem-fpkm-collapsed.polya.rds --output gsva_scores_stranded.tsv`


+ **Eight** files named as `results/gsva_<tukey/anova>_<stranded/polya>_<disease_type_new/short_histology>.tsv` represent results from modeling
	+ Files created with: `Rscript --vanilla 02-model-gsea.R`
	+ Assumes `results/gsva_scores_stranded.tsv` and `results/gsva_scores_polya.tsv` exist
 