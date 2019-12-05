#!/bin/bash
# JA Shapiro
# CCDL for ALSF 2019
set -e
set -o pipefail

COSMIC_USER=${OPENPBTA_COSMIC_USER:-""}
COSMIC_PASS${OPENPBTA_COSMIC_PASS:-""}


# The sqlite database made from the callers will be called:
dbfile=scratch/vus_db.sqlite
consensus=data/pbta-snv-lancet.vep.maf.gz
metadata=data/pbta-histologies.tsv

python3 analyses/recurrent-VUS/scripts/01-setup_db.py \
  --db-file $dbfile \
  --consensus-file $consensus \
  --meta-file $metadata
  --cosmic-user $COSMIC_USER
  --cosmic-pass $COSMIC_PASS