This directory contains CSV files that represent data inputted to plots shown in the manuscript.
Data for the following figures is included.

| Figure | Brief Description of Figure |
|--------|-------------|
|   1B     |    Barplot showing the number of biospecimens per phase of therapy       |
|   2A     |    Onocoprint summarizing low-grade glioma alterations across samples for top 20 mutated genes          |
|   2B     |    Onocoprint summarizing embryonal tumor alterations across samples for top 20 mutated genes         |
|   2C     |    Onocoprint summarizing high-grade glioma alterations across samples for top 20 mutated genes         |
|   2D     |    Onocoprint summarizing other CNS tumor alterations across samples for top 20 mutated genes         |
|   3A     |    Barplot of occurrence and co-occurrence of nonsynonymous mutations for the 50 most commonly mutated genes across all tumor types         |
|   3B     |  Heatmap of co-occurrence and mutual exclusivity of nonsynonymous mutations between genes           |
|   3C     |  Scatterplot showing correlation of SV and CNV breaks            |
|   3D     |  Barplot showing chromothripsis frequency for all cancer groups with N >= 3 tumors           |
|   3E     |  Sina plots of RefSig signature weights across cancer groups         |
|   4A     |  ROC curve for TP53 classifier         |
|   4B     |  Violin/strip plots of TP53 scores across TP53 alteration status           |
|   4C     |   Violin/strip plots of TP53 expression [log(FPKM)] across TP53 alteration status            |
|   4D     |  Box/strip plots of TP53 scores and telomerase scores across cancer groups          |
|   4E    |    Heatmap of RegSg signatures in patients with a hypermutator phenotype       |
|   4F     |  Forest plot from survival analysis exploring TP53 and telomerase score effects          |
|   4G     |  Forest plot from survival analysis exploring HGG molecular subtype effects          |
|   4H     |  Kaplan-Meier cure of HGG tumors by molecular subtype           |
|   5A     |  UMAP of tumors colored by broad histology           |
|   5B     |  Heatmap of GSVA scores across cancer groups           |
|   5C     |  Box/strip plots of quanTIseq scores across cancer groups           |
|   5D     |  Forest plot from survival analysis exploring CD274 expression and immune cell proportion effects           |
|   5E     |  Box/strip plot of CD274 expression across medulloblastoma subtypes          |
|   S2A     |             |
|   S2B     |             |
|   S2C     |             |
|   S3D     |             |
|   S2E     |             |
|   S2F     |             |
|   S2G     |             |
|   S2H     |             |
|   S2I     |             |
|   S3A     |             |
|   S3B     |             |
|   S3C     |             |
|   S3D     |             |
|   S3E     |             |
|   S4A     |             |
|   S4B     |             |
|   S5A     |             |
|   S5B     |             |
|   S5C     |             |
|   S6A     |             |
|   S6B     |             |
|   S6C     |             |
|   SD6     |             |
|   S6E     |             |
|   S6F     |             |
|   S7A     |             |
|   S7B     |             |
|   S7C     |             |
|   S7D     |             |
|   S7E     |             |
|   S7F     |             |
|   S7G     |             |
|   S7H     |             |
|   S7I     |             |






|--------|--------|------------------|
| [`fig1-sample-distribution.R`](fig1-sample-distribution.R) | Figure 1B | N/A |
| [`fig2_figS3-oncoprint-landscape.R`](fig2-oncoprint-landscape.R) | Figure 2 & Figure S3B | [`oncoprint-landscape`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/oncoprint-landscape) |
| [`fig3-chromothripsis-barplot.R`](fig3-chromothripsis-barplot.R) | Figure 3D | [`chromothripsis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromothripsis) |
| [`fig3_figS4-mutational-signatures-panels.R`](fig3_figS4-mutational-signatures-panels.R) | Figure 3E and Figure S4 | [`mutational-signatures`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures) |
| [`fig4-tp53-telomerase-panels.R`](fig4-tp53-telomerase-panels.R) | Figure 4 panels A, B, C, D, F | [`tp53_nf1_score`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tp53_nf1_score) <br> [`telomerase-activity-prediction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction/) <br> [`survival-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/survival-analysis)  |
| [`fig4-heatmap.R`](fig4-heatmap.R) | Figure 4E | [`mutational-signatures`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/mutational-signatures)  |
| [`fig4-hgg-kaplan-meier.R`](fig4-hgg-kaplan-meier.R) | Figure 4H | [`survival-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/survival-analysis)  |
| [`fig4-hgg-subtype-forest-plot.R`](fig4-hgg-subtype-forest-plot.R) | Figure 4G | [`survival-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/survival-analysis)  |
| [`fig5-panels-gsva-umap.R`](fig5-panels-gsva-umap.R) | Figure 5 panels A, B | [`transcriptomic-dimension-reduction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/transcriptomic-dimension-reduction) <br> [`collapse-rnaseq`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/collapse-rnaseq) <br> [`gene-set-enrichment-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/gene-set-enrichment-analysis) |
| [`fig5-forest-plot.R`](fig5-forest-plot.R) | Figure 5D | [`survival-analysis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/survival-analysis) <br> [`immune-deconv`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/immune-deconv)
| [`fig5_figS6-immune-deconv-panels.R`](fig5_figS6-immune-deconv-panels.R) |  Figure 5 panels C, E and Figure S6 panels E, F | [`immune-deconv`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/immune-deconv)  |
| [`figS2-snv-callers-panels.R`](figS2-snv-callers-panels.R) | Figure S2 panels A-G | [`snv-callers`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/snv-callers) (:warning: _Requires 32 GB RAM_)
| [`figS2-tmb-compare-panels.R`](figS2-tmb-compare-panels.R) | Figure S2 panels H, I | [`tmb-compare`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tmb-compare)
| [`figS3-panels-cn-chromothripsis.R`](figS3-panels-cn-chromothripsis.R) | Figure S3 panels C, D, E | [`cnv-chrom-plot`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/cnv-chrom-plot) <br> [`chromothripsis`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromothripsis) <br> [`copy_number_consensus_call`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/copy_number_consensus_call)
| [`figS5-all-panels.R`](figS5-all-panels.R) | Figure S5 A-C | [`tp53_nf1_score`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tp53_nf1_score) <br> [`telomerase-activity-prediction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction/)
| [`figS6-subtype-umap-panels.R`](figS6-subtype-umap-panels.R) | Figure S6 panels A-D | [`transcriptomic-dimension-reduction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/transcriptomic-dimension-reduction)
| [`figS7-seqcenter-barplot.R`](figS7-seqcenter-barplot.R) | Figure S7 panel A | None |
| [`figS7-UMAP-libraries.R`](figS7-UMAP-libraries.R) | Figure S7 panel B | None |
| [`figS7-tp53-telomerase-tumor-purity-threshold.R`](figS7-tp53-telomerase-tumor-purity-threshold.R) | Figure S7 panel G | [`transcriptomic-dimension-reduction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/transcriptomic-dimension-reduction) <br> [`telomerase-activity-prediction`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction/) <br> [`tumor-purity-exploration`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/tumor-purity-exploration/)
