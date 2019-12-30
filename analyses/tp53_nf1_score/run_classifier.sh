#!/bin/bash

# K S Gaonkar

collapsed_stranded="pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
collapsed_polya="pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"

# Run classifier for stranded and polya
python3 analyses/tp53_nf1_score/01-apply-classifier.py -f $collapsed_stranded
python3 analyses/tp53_nf1_score/01-apply-classifier.py -f $collapsed_polya


