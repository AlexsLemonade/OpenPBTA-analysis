#!/bin/bash

set -e
set -o pipefail

python3 analyses/snv-callers/scripts/01-setup_db.py \
  --db-file scratch/snv_db.sqlite \
  --strelka-file data/pbta-snv-strelka2.vep.maf.gz \
  --mutect-file data/pbta-snv-mutect2.vep.maf.gz \
  --lancet-file data/pbta-snv-lancet.vep.maf.gz \
  --meta-file data/pbta-histologies.tsv
