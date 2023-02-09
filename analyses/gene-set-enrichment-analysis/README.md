# Gene Set Enrichment Analysis

Written by Stephanie J. Spielman to supercede previous analyses in [`ssgsea-hallmark`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/ssgsea-hallmark).

Primary goals include:

1. Score hallmark pathways based on expression data using GSVA analysis, using a strategy that produces Gaussian-distributed scores.
2. Analyze scores for highly significant differences among tumor classifications

## Usage:

Note that running this analyis on the full dataset requires > 16GB of memory.
Run the bash script of this analysis module:

using OPENPBTA_BASE_SUBTYPING=1 to run this module using the pbta-histologies-base.tsv from data folder while running molecular-subtyping modules for release.
```sh
OPENPBTA_BASE_SUBTYPING=1 analyses/gene-set-enrichment-analysis/run-gsea.sh
```

OR by default uses pbta-histologies.tsv from data folder
```sh
bash analyses/gene-set-enrichment-analysis/run-gsea.sh
```

*This command above assumes you are in the top directory, OpenPBTA-analysis*

## Folder Content

+ `01-conduct-gsea-analysis.R` performs the GSVA analysis using RSEM FPKM expression data for both stranded and polyA data.
Results are saved in `results/` TSV files when run via `run-gsea.sh`.

+ `02-model-gsea.Rmd` performs ANOVA and Tukey tests on GSVA scores to evaluate, for each hallmark pathway, differences in GSVA across groups (e.g. short histology or disease type).

+ `03-assess-gsea-at-threshold.Rmd` compares certain key GSEA results generated from the full stranded data compared to stranded data that passes a tumor purity threshold.

+ `results/gsva_scores_stranded.tsv` represents GSVA scores calculated from `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` with `Rscript --vanilla 01-conduct-gsea-analysis.R`

+ `results/gsva_scores_polya.tsv` represents GSVA scores calculated from `pbta-gene-expression-rsem-fpkm-collapsed.polya.rds` with with `Rscript --vanilla 01-conduct-gsea-analysis.R`

+ `results/gsva_scores_stranded_threholded.tsv` represents GSVA scores calculated from `pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds` with `Rscript --vanilla 01-conduct-gsea-analysis.R --apply_tumor_purity_threshold`

+ **Eight** files named as `results/gsva_<tukey/anova>_<stranded/polya>_<harmonized_diagnosis/cancer_group>.tsv` represent results from modeling
	+ Files created with: `Rscript --vanilla 02-model-gsea.R`
	+ Assumes `results/gsva_scores_stranded.tsv` and `results/gsva_scores_polya.tsv` exist

