#!/bin/bash

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples

# The script takes one environment variable, `OPENPBTA_BASE_SUBTYPING`, if value is 1 then
# uses histologies-base.tsv for subtyping if value is 0 runs all modules with histologies.tsv(Default)

set -e
set -o pipefail

RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
analysis_dir="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$analysis_dir" || exit


# we want to skip the poly-A steps in CI
# if POLYA=1, poly-A steps will be run
POLYA=${OPENPBTA_POLYAPLOT:-1}

data_dir="../../data"
scratch_dir="../../scratch"
# cds gencode bed file  
cds_file="${scratch_dir}/gencode.v27.primary_assembly.annotation.bed"
snvconsensus_file="${data_dir}/snv-consensus-plus-hotspots.maf.tsv.gz"
cnvconsensus_file="${data_dir}/consensus_wgs_plus_cnvkit_wxs.tsv.gz"

if [[ RUN_FOR_SUBTYPING == "0" ]]
then
   histology_file="../../data/histologies.tsv" 
else 
   histology_file="../../data/histologies-base.tsv"  
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

if [[ RUN_FOR_SUBTYPING == "0" ]]
then
   # expression files for prediction
   collapsed_rna="${data_dir}/gene-expression-rsem-tpm-collapsed.rds"
else
   # expression files for prediction
   collapsed_rna="../collapse-rnaseq/results/gene-expression-rsem-tpm-collapsed.rds"
fi


# Run classifier and ROC plotting for RNA data
python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_rna} -t ${histology_file}

# check correlation expression and scores
Rscript -e "rmarkdown::render('${analysis_dir}/02-qc-rna_expression_score.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# subset cnv where tp53 is lost
Rscript -e "rmarkdown::render('${analysis_dir}/03-tp53-cnv-loss-domain.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# subset SV where tp53 is lost
Rscript -e "rmarkdown::render('${analysis_dir}/04-tp53-sv-loss.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# gather TP53 altered status
Rscript -e "rmarkdown::render('${analysis_dir}/05-tp53-altered-annotation.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"

# evaluate classifer scores for stranded data
python3 ${analysis_dir}/06-evaluate-classifier.py -s ${analysis_dir}/results/tp53_altered_status.tsv -f ${analysis_dir}/results/gene-expression-rsem-tpm-collapsed_classifier_scores.tsv -c ${histology_file} 

# Skip poly-A steps in CI
if [ "$POLYA" -gt "0" ]; then
  python3 ${analysis_dir}/06-evaluate-classifier.py -s ${analysis_dir}/results/tp53_altered_status.tsv -f ${analysis_dir}/results/gene-expression-rsem-tpm-collapsed_classifier_scores.tsv -c ${histology_file}
fi

