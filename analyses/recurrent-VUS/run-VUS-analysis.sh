#!/bin/bash
# JA Shapiro
# CCDL for ALSF 2019
set -e
set -o pipefail

# The sqlite database made from the callers will be called:
dbfile=scratch/vus_db.sqlite
consensus=data/pbta-snv-lancet.vep.maf.gz
metadata=data/pbta-histologies.tsv

python3 analyses/recurrent-VUS/scripts/01-setup_db.py \
  --db-file $dbfile \
  --consensus-file $consensus
  --meta-file $metadata