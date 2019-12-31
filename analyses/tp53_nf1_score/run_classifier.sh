#!/bin/bash

set -e
set -o pipefail

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples


collapsed_stranded="pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
collapsed_polya="pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"

# Run classifier for stranded and polya
python3 analyses/tp53_nf1_score/01-apply-classifier.py -f $collapsed_stranded
python3 analyses/tp53_nf1_score/01-apply-classifier.py -f $collapsed_polya


