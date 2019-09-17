# Bethell and Taroni for CCDL 2019
# Run dimension reduction for all subsets of gene expression for both methods:
# RSEM and kallisto
# 
# Usage: bash 02-dimension-reduction.sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

SEED=2019
PERPLEXITY=10
NEIGHBORS=15
COUNT_THRESHOLD=100
METADATA="../../data/pbta-histologies.tsv"
OUTPUT="results"

#### RSEM ----------------------------------------------------------------------

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-rsem.fpkm.rds \
  --metadata ${METADATA} \
  --filename_lead rsem_all \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../scratch/pbta-gene-expression-rsem_stranded.fpkm.rds \
  --metadata ${METADATA} \
  --filename_lead rsem_stranded \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../scratch/pbta-gene-expression-rsem_polyA.fpkm.rds \
  --metadata ${METADATA} \
  --filename_lead rsem_polyA \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}
  
#### kallisto ------------------------------------------------------------------

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-kallisto.rds \
  --metadata ${METADATA} \
  --filename_lead kallisto_all \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../scratch/pbta-gene-expression-kallisto_stranded.rds \
  --metadata ${METADATA} \
  --filename_lead kallisto_stranded \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../scratch/pbta-gene-expression-kallisto_polyA.rds \
  --metadata ${METADATA} \
  --filename_lead kallisto_polyA \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}
