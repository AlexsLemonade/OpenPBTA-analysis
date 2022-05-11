#!/bin/sh

set -e
set -o pipefail

# This script runs RNA-seq summary modules 
# when new RNA samples are added to OpenPBTA
# Subtyping modules then need to be re-run for 
# these new samples when needed for data release

## Step 1. Generate summary files needed for subtyping that are not included in data download

echo "Run GSVA"
OPENPBTA_BASE_SUBTYPING=1 bash ../analyses/gene-set-enrichment-analysis/run-gsea.sh

echo "TP53 altered score"
OPENPBTA_BASE_SUBTYPING=1 bash ../analyses/tp53_nf1_score/run_classifier.sh

echo "Chromosomal instability"
OPENPBTA_BASE_SUBTYPING=1 bash ../analyses/chromosomal-instability/run_breakpoint_analysis.sh

## Step 2. Run subtyping modules

echo "Run MB subtyping"
bash ../analyses/molecular-subtyping-MB/run-molecular-subtyping-mb.sh

echo "Run CRANIO subtyping"
bash ../analyses/molecular-subtyping-CRANIO/run-molecular-subtyping-cranio.sh

echo "Run EPN subtyping"
bash ../analyses/molecular-subtyping-EPN/run-molecular-subtyping-EPN.sh

echo "Run Embryonal subtyping"
bash ../analyses/molecular-subtyping-embryonal/run-embryonal-subtyping.sh

echo "Run chordoma subtyping"
bash ../analyses/molecular-subtyping-chordoma/run-molecular-subtyping-chordoma.sh

echo "Run EWS subtyping"
bash ../analyses/molecular-subtyping-EWS/run_subtyping.sh

echo "Run neurocytoma subtyping"
bash ../analyses/molecular-subtyping-neurocytoma/run_subtyping.sh

echo "Run HGG subtyping"
bash ../analyses/molecular-subtyping-HGG/run-molecular-subtyping-HGG.sh

echo "Run LGAT subtyping"
bash ../analyses/molecular-subtyping-LGAT/run_subtyping.sh

echo "Run compile subtyping"
bash ../analyses/molecular-subtyping-pathology/run-subtyping-aggregation.sh

echo "Run integrate subtyping"
bash ../analyses/molecular-subtyping-integrate/run-subtyping-integrate.sh
