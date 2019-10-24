#!/bin/bash

# K S Gaonkar 

# Run fusion_filtering

set -e
set -o pipefail

# Run Fusion standardization
Rscript analyses/fusion_filtering/01-fusion-standardization.R -f "data/pbta-fusion-arriba.tsv.gz" -c "arriba" -o "scratch/arriba.tsv"
Rscript analyses/fusion_filtering/01-fusion-standardization.R -f "data/pbta-fusion-starfusion.tsv.gz" -c "starfusion" -o "scratch/starfusion.tsv"

# Run Fusion general filtering for arriba
Rscript analyses/fusion_filtering/02-fusion-filtering.R -S scratch/arriba.tsv -e data/pbta-gene-expression-rsem-fpkm.polya.rds -r -a "GTEx|HGNC_GENEFAM|DGD_PARALOGS|Normal|BodyMap|ConjoinG" -j 1 -s 10 -i "in-frame|frameshift|other" -R analyses/fusion_filtering/references/ -o scratch/arribaPolya -t 1

Rscript analyses/fusion_filtering/02-fusion-filtering.R -S scratch/arriba.tsv -e data/pbta-gene-expression-rsem-fpkm.stranded.rds -r -a "GTEx|HGNC_GENEFAM|DGD_PARALOGS|Normal|BodyMap|ConjoinG" -j 1 -s 10 -i "in-frame|frameshift|other" -R analyses/fusion_filtering/references/ -o scratch/arribaStranded -t 1

# Run Fusion general filtering for atarfusion
Rscript analyses/fusion_filtering/02-fusion-filtering.R -S scratch/starfusion.tsv -e data/pbta-gene-expression-rsem-fpkm.polya.rds -r -a "GTEx|HGNC_GENEFAM|DGD_PARALOGS|Normal|BodyMap|ConjoinG" -j 1 -s 10 -i "in-frame|frameshift|other" -R analyses/fusion_filtering/references/ -o scratch/starfusionPolya -t 1

Rscript analyses/fusion_filtering/02-fusion-filtering.R -S scratch/starfusion.tsv -e data/pbta-gene-expression-rsem-fpkm.stranded.rds -r -a "GTEx|HGNC_GENEFAM|DGD_PARALOGS|Normal|BodyMap|ConjoinG" -j 1 -s 10 -i "in-frame|frameshift|other" -R analyses/fusion_filtering/references/ -o scratch/starfusionStranded -t 1

# Fusion Annotation
Rscript analyses/fusion_filtering/03-Calc-zscore-annotate.R -S scratch/arribaPolya_QC_expression_filtered_annotated.RDS -e data/pbta-gene-expression-rsem-fpkm.polya.rds -g analyses/fusion_filtering/references/Brain_FPKM_hg38_matrix.txt.zip -o scratch/arribaPolya_QC_expression 

Rscript analyses/fusion_filtering/03-Calc-zscore-annotate.R -S scratch/arribaStranded_QC_expression_filtered_annotated.RDS -e data/pbta-gene-expression-rsem-fpkm.polya.rds -g analyses/fusion_filtering/references/Brain_FPKM_hg38_matrix.txt.zip -o scratch/arribaStranded_QC_expression

Rscript analyses/fusion_filtering/03-Calc-zscore-annotate.R -S scratch/starfusionPolya_QC_expression_filtered_annotated.RDS -e data/pbta-gene-expression-rsem-fpkm.polya.rds -g analyses/fusion_filtering/references/Brain_FPKM_hg38_matrix.txt.zip -o scratch/starfusionPolya_QC_expression

Rscript analyses/fusion_filtering/03-Calc-zscore-annotate.R -S scratch/starfusionStranded_QC_expression_filtered_annotated.RDS -e data/pbta-gene-expression-rsem-fpkm.polya.rds -g analyses/fusion_filtering/references/Brain_FPKM_hg38_matrix.txt.zip -o scratch/starfusionStranded_QC_expression
