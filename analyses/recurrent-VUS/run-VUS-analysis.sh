#!/bin/bash
# JA Shapiro
# CCDL for ALSF 2019
set -e
set -o pipefail

# The sqlite database made from the callers will be called:
dbfile=scratch/vus_db.sqlite
consensus=data/snv-consensus_11122019/consensus_mutation.maf.tsv
metadata=data/pbta-histologies.tsv

python3 analyses/recurrent-VUS/scripts/01-setup_db.py \
  --db-file $dbfile \
  --conensus-file $consensus
  --meta-file $metadata