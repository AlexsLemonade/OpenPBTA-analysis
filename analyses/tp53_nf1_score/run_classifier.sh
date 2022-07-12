#!/bin/bash

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples

# The script can tak three environment variables:
# `OPENPBTA_BASE_SUBTYPING`: 
#     if value is 1, then uses pbta-histologies-base.tsv for subtyping. 
#     if value is 0 (DEFAULT), runs module with pbta-histologies.tsv
# `OPENPBTA_TP53_FIGURES`: 
#     if value is 1, uses expression data in the data release for prediction
#     if value is 0 (DEFAULT), uses expression data in the `collapse-rnaseq` module for prediction
# `OPENPBTA_POLYAPLOT`: 
#     if value is 1 (DEFAULT), runs the POLYA steps
#     if value is 0, skips the POLYA steps

set -e
set -o pipefail

# If running for subtyping, use the base histologies file by setting this to 1
RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

# If we're running for figure generation, set to 1
RUN_FOR_FIGURES=${OPENPBTA_TP53_FIGURES:-0}

# we want to skip the poly-A steps in CI
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

# cds gencode bed file
cds_file="${scratch_dir}/gencode.v27.primary_assembly.annotation.bed"
snvconsensus_file="${data_dir}/pbta-snv-consensus-mutation.maf.tsv.gz"
cnvconsensus_file="${data_dir}/consensus_seg_annotated_cn_autosomes.tsv.gz"

if [[ "$RUN_FOR_SUBTYPING" -eq "0" ]]
then
   histology_file="../../data/pbta-histologies.tsv"
else
   histology_file="../../data/pbta-histologies-base.tsv"
fi


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

# If running for the purpose of figure generation (RUN_FOR_FIGURES will be 1),
# use the data in the analysis directory that has been freshly collapsed
if [[ "$RUN_FOR_FIGURES" -eq "0" ]]
then
  # expression files for prediction
  collapsed_stranded="${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
  collapsed_polya="${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"
else
  # expression files for prediction
  collapsed_stranded="../collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
  collapsed_polya="../collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"
fi

# Skip poly-A steps in CI
if [ "$POLYA" -gt "0" ]; then
  python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_polya}
fi

# Run classifier and ROC plotting for stranded data
python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_stranded}

# check correlation expression and scores
Rscript -e "rmarkdown::render('${analysis_dir}/02-qc-rna_expression_score.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# subset cnv where tp53 is lost
Rscript -e "rmarkdown::render('${analysis_dir}/03-tp53-cnv-loss-domain.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# subset SV where tp53 is lost
Rscript -e "rmarkdown::render('${analysis_dir}/04-tp53-sv-loss.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# gather TP53 altered status
Rscript -e "rmarkdown::render('${analysis_dir}/05-tp53-altered-annotation.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# evaluate classifer scores for stranded data
python3 ${analysis_dir}/06-evaluate-classifier.py -s ${analysis_dir}/results/tp53_altered_status.tsv -f ${analysis_dir}/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded_classifier_scores.tsv -c ${histology_file} -o stranded

# Skip poly-A steps in CI
if [ "$POLYA" -gt "0" ]; then
  python3 ${analysis_dir}/06-evaluate-classifier.py -s ${analysis_dir}/results/tp53_altered_status.tsv -f ${analysis_dir}/results/pbta-gene-expression-rsem-fpkm-collapsed.polya_classifier_scores.tsv -c ${histology_file} -o polya
fi

# plot ROC curves for poly-A and stranded data
Rscript 07-plot-roc.R

# Skip plotting steps for subtyping
if [[ "$RUN_FOR_SUBTYPING" -eq "0" ]]; then
# create violin plots of TP53 scores across molecular subtypes per broad histology
  Rscript 08-compare-molecularsubtypes-tp53scores.R

# create boxplots by broad histology, cancer group
  Rscript 09-compare-histologies.R
fi
