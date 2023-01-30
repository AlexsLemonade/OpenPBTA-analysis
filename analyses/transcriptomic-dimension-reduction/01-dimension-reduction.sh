# Bethell and Taroni for CCDL 2019
# Run dimension reduction for all subsets of gene expression for both methods:
# RSEM and kallisto
#
# Usage: bash 01-dimension-reduction.sh

# Takes one environment variable, `BASE_SUBTYPING`, if value is 1 then
# uses pbta-histologies-base.tsv for subtyping if value is 0 runs all modules with pbta-histologies.tsv(Default)


# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

SEED=2019
PERPLEXITY=10
NEIGHBORS=15
COUNT_THRESHOLD=100
if [ "$RUN_FOR_SUBTYPING" == 0 ]; then
METADATA="../../data/pbta-histologies.tsv"
else
METADATA="../../data/pbta-histologies-base.tsv" 
fi

OUTPUT="results"

#### RSEM ----------------------------------------------------------------------

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-rsem-fpkm.polya.rds \
  --metadata ${METADATA} \
  --filename_lead rsem_polyA_none \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-rsem-fpkm.polya.rds \
  --metadata ${METADATA} \
  --filename_lead rsem_polyA_log \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD} \
  --log2_transform

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-rsem-fpkm.stranded.rds \
  --metadata ${METADATA} \
  --filename_lead rsem_stranded_none \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-rsem-fpkm.stranded.rds \
  --metadata ${METADATA} \
  --filename_lead rsem_stranded_log \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD} \
  --log2_transform
  
  

#### RSEM log2 stranded, excluding mitochondrial genes and skipping t-sne ------------------------------------
# Note that these results are only used for exploration in `05-seq_center-mitochondrial-genes.Rmd`
Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-rsem-fpkm.stranded.rds \
  --metadata ${METADATA} \
  --remove_mito_genes \
  --filename_lead rsem_stranded_no-mito_log \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD} \
  --log2_transform \
  --skip_tsne
  
  
#### kallisto ------------------------------------------------------------------

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-kallisto.polya.rds \
  --metadata ${METADATA} \
  --filename_lead kallisto_polyA_none \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-kallisto.polya.rds \
  --metadata ${METADATA} \
  --filename_lead kallisto_polyA_log \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD} \
  --log2_transform

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-kallisto.stranded.rds \
  --metadata ${METADATA} \
  --filename_lead kallisto_stranded_none \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD}

Rscript --vanilla scripts/run-dimension-reduction.R \
  --expression ../../data/pbta-gene-expression-kallisto.stranded.rds \
  --metadata ${METADATA} \
  --filename_lead kallisto_stranded_log \
  --output_directory ${OUTPUT} \
  --seed ${SEED} \
  --perplexity ${PERPLEXITY} \
  --neighbors ${NEIGHBORS} \
  --low_count_threshold ${COUNT_THRESHOLD} \
  --log2_transform



