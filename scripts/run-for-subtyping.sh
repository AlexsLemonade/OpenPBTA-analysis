#!/bin/sh

set -e
set -o pipefail

# This script runs RNA-seq summary modules 
# when new RNA samples are added to OpenPBTA
# Subtyping modules then need to be re-run for 
# these new samples when needed for data release

BASE_SUBTYPING=1

## Step 1. Generate summary files needed for subtyping

echo "Create collapse rsem files"
BASE_SUBTYPING=1 ../analyses/collapse-rnaseq/run-collapse-rnaseq.sh

echo "Create independent sample list"
BASE_SUBTYPING=1 ../analyses/independent-samples/run-independent-samples.sh 

echo "Create fusion filtered list" 
BASE_SUBTYPING=1 ../analyses/fusion_filtering/run_fusion_merged.sh

echo "Run fusion summary for subtypes"
BASE_SUBTYPING=1 ../analyses/fusion-summary/run-new-analysis.sh

echo "Run tsne"
BASE_SUBTYPING=1 ../analyses/transcriptomic-dimension-reduction/dimension-reduction-plots.sh

echo "Run gsea"
BASE_SUBTYPING=1 ../analyses/gene-set-enrichment-analysis/run-gsea.sh



## Step 2. Run subtyping modules

echo "Run MB subtyping"
bash ../analyses/molecular-subtyping-MB/run-molecular-subtyping-mb.sh

echo "Run CRANIO subtyping"
bash ../analyses/molecular-subtyping-CRANIO/run-molecular-subtyping-cranio.sh

echo "Run EPN subtyping"
bash ../analyses/molecular-subtyping-EPN/run-molecular-subtyping-EPN.sh

echo "Run Embryonal subtyping"
bash ../analyses/molecular-subtyping-embryonal/run-embryonal-subtyping.sh

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

