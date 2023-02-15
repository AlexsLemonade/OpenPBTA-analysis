#!/bin/bash

# SJ Spielman 2023 for CCDL

# This script runs the relevant aspects of the classifier pipeline (see [`run_classifier.sh`](./run_classifier.sh)) on tumors which have passed the tumor purity threshold
# Note that certain notebooks/scripts from the module are not run here since they are not strictly for this re-analysis:
# - 02-qc-rna_expression_score.Rmd
# - 07-plot-roc.R
# - 08-compare-molecularsubtypes-tp53scores.R
# - 09-compare-histologies.R
# The scripts/notebooks we do run here all are needed as they produce outputs to be consumed by the additional notebook, `10-tp53-tumor-purity-threshold.Rmd`,
#  which performs certain manuscript-level assessments on the thresholded analysis (including ROC).
# Note further that this script should be run _after_ subtyping, not during. PolyA-specific scripts do not need to be run here.

set -euo pipefail

# We want to skip the poly-A steps in CI
# if POLYA=1, poly-A steps will be run
POLYA=${OPENPBTA_POLYAPLOT:-1}


# This script should always run as if it were being called from
# the directory it lives in.
analysis_dir="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$analysis_dir" || exit

data_dir="../../data"
scratch_dir="../../scratch/tp53-classifier"
# Make sure scratch directory exists
mkdir -p $scratch_dir

# Output directory for tumor purity threshold results and HTML notebooks
output_dir="${analysis_dir}/results/tumor-purity-threshold"
mkdir -p $output_dir

# cds gencode bed file
cds_file="${scratch_dir}/gencode.v27.primary_assembly.annotation.bed"
snvconsensus_file="${data_dir}/pbta-snv-consensus-mutation.maf.tsv.gz"
cnvconsensus_file="${data_dir}/consensus_seg_annotated_cn_autosomes.tsv.gz"

histology_file="../../data/pbta-histologies.tsv"


# Step 1: Prep the SNV consensus data for evaluation downstream
# This only needs to be run if the output file `TP53_NF1_snv_alteration.tsv` does not already exist.
if [ ! -f  ${analysis_dir}/results/TP53_NF1_snv_alteration.tsv ]; then
  # Convert GTF to BED file
  # Here we are only extracting lines with as a CDS i.e. are coded in protein
  gunzip -c ${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz \
    | awk '$3 ~ /CDS/' \
    | convert2bed --do-not-sort --input=gtf - \
    > $cds_file

  # Prep the SNV consensus data for evaluation downstream
  Rscript --vanilla ${analysis_dir}/00-tp53-nf1-alterations.R \
    --snvConsensus ${snvconsensus_file} \
    --cnvConsensus ${cnvconsensus_file} \
    --histologyFile ${histology_file} \
    --outputFolder ${analysis_dir}/results \
    --gencode ${cds_file}
fi

# expression files for prediction from data release
collapsed_stranded="${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"

# Run classifier and ROC plotting for stranded data, with tumor purity threshold turned on (`-t` flag)
IDS_FILE="${analysis_dir}/../tumor-purity-exploration/results/thresholded_rna_stranded_same-extraction.tsv"  # file for filtering IDs
python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_stranded} -t --ids_file ${IDS_FILE}

# Run notebook to produce the `loss_overlap_domains_tp53` file consumed by 05-tp53-altered-annotation.Rmd
OUTPUT_HTML_03=${output_dir}/03-tp53-cnv-loss-domain_tumor-purity-threshold.nb.html # output HTML file name
Rscript -e "rmarkdown::render('${analysis_dir}/03-tp53-cnv-loss-domain.Rmd', output_file = '${OUTPUT_HTML_03}', params=list(tumor_purity_threshold = 1))"

# Run notebook to produce the `fusion_bk_tp53_loss` and `sv_overlap_tp53` files consumed by 05-tp53-altered-annotation.Rmd
OUTPUT_HTML_04=${output_dir}/04-tp53-sv-loss_tumor-purity-threshold.nb.html # output HTML file name
Rscript -e "rmarkdown::render('${analysis_dir}/04-tp53-sv-loss.Rmd', output_file = '${OUTPUT_HTML_04}', params=list(tumor_purity_threshold = 1))"

# Run notebook to produce `tp53_altered_status` file consumed by 06-evaluate-classifier.py
OUTPUT_HTML_05=${output_dir}/05-tp53-altered-annotation_tumor-purity-threshold.nb.html # output HTML file name
Rscript -e "rmarkdown::render('${analysis_dir}/05-tp53-altered-annotation.Rmd', output_file = '${OUTPUT_HTML_05}', params=list(tumor_purity_threshold = 1))"

# Produce files needed for ROC
python3 ${analysis_dir}/06-evaluate-classifier.py \
  -s ${analysis_dir}/results/tp53_altered_status_tumor-purity-threshold.tsv \
  -f ${analysis_dir}/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded_classifier_scores_tumor-purity-threshold.tsv \
  -c ${histology_file} \
  -o stranded_tumor-purity-threshold \
  -r ${output_dir} \
  -p ${output_dir}
