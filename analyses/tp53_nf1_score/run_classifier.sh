#!/bin/bash

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples

set -e
set -o pipefail

# we want to skip the poly-A steps in CI
# if POLYA=1, poly-A steps will be run
POLYA=${OPENPBTA_POLYAPLOT:-1}

data_dir="data"
scratch_dir="scratch"
# cds gencode bed file  
cds_file="${scratch_dir}/gencode.v27.primary_assembly.annotation.bed"
snvconsensus_file="${data_dir}/pbta-snv-consensus-mutation.maf.tsv.gz"
cnvconsensus_file="${data_dir}/consensus_seg_annotated_cn_autosomes.tsv.gz"
histology_file="${data_dir}/pbta-histologies.tsv"
analysis_dir="analyses/tp53_nf1_score"

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

# expression files for prediction
collapsed_stranded="pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
collapsed_polya="pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"


# Skip poly-A steps in CI
if [ "$POLYA" -gt "0" ]; then
  python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_polya}
fi

# Run classifier and ROC plotting for stranded data
python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_stranded}

# check correlation expression and scores
Rscript -e "rmarkdown::render('${analysis_dir}/02-qc-rna_expression_score.Rmd')"

# subset cnv where tp53 is lost
Rscript -e "rmarkdown::render('${analysis_dir}/03-tp53-cnv-loss-domain.Rmd')"

# subset SV where tp53 is lost
Rscript -e "rmarkdown::render('${analysis_dir}/04-tp53-sv-loss.Rmd')"

# gather TP53 altered status
Rscript -e "rmarkdown::render('${analysis_dir}/05-tp53-altered-annotation.Rmd')"

# evaluate classifer scores for stranded data
python3 ${analysis_dir}/06-evaluate-classifier.py -s ${analysis_dir}/results/tp53_altered_status.tsv -f ${analysis_dir}/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded_classifier_scores.tsv -c ${data_dir}/pbta-histologies.tsv -o stranded

# Skip poly-A steps in CI
if [ "$POLYA" -gt "0" ]; then
  python3 ${analysis_dir}/06-evaluate-classifier.py -s ${analysis_dir}/results/tp53_altered_status.tsv -f ${analysis_dir}/results/pbta-gene-expression-rsem-fpkm-collapsed.polya_classifier_scores.tsv -c ${data_dir}/pbta-histologies.tsv -o polya
fi

