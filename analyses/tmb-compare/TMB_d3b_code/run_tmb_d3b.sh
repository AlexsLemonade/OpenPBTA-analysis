#!/bin/bash

# fail on error
set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

BASE_DIR=../../../

python3 code/01_calculate_tmb.py \
  -i $BASE_DIR/data/pbta-snv-consensus-mutation.maf.tsv.gz \
  -m $BASE_DIR/data/pbta-histologies.tsv \
  -c config_files/calculate_tmb.cfg.json \
  -w config_files/target_cfg.targetcombos.txt \
  -d inputs/target_files \
  -o outputs/pbta-snv-consensus-TMB_intarget.txt

python3 code/02_cumulative_freq_TMBplot.py  \
  -t outputs/pbta-snv-consensus-TMB_intarget.txt \
  -o outputs/pbta-snv-consensus \
  -s 10
