#!/bin/bash

set -e
set -o pipefail

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples

collapsed_stranded="pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
collapsed_polya="pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"

#cds gencode bed file  
exon_file="scratch/gencode.v27.primary_assembly.annotation.bed"

# Run classifier for stranded and polya
python3 analyses/tp53_nf1_score/01-apply-classifier.py -f $collapsed_stranded
python3 analyses/tp53_nf1_score/01-apply-classifier.py -f $collapsed_polya

# Convert GTF to BED file for use in bedtools
# Here we are only extracting lines with as a CDS i.e. are coded in protein
gunzip -c data/gencode.v27.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  > $exon_file
