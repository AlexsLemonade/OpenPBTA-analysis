#!/bin/bash

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples

set -e
set -o pipefail

# we want to skip the poly-A ROC plot in CI
POLYA_PLOT=${OPENPBTA_POLYAPLOT:-1}

data_dir="data"
scratch_dir="scratch"
# cds gencode bed file  
cds_file="${scratch_dir}/gencode.v27.primary_assembly.annotation.bed"
consensus_file="${data_dir}/pbta-snv-consensus-mutation.maf.tsv.gz"
clinical_file="${data_dir}/pbta-histologies.tsv"
analysis_dir="analyses/tp53_nf1_score"

# Convert GTF to BED file
# Here we are only extracting lines with as a CDS i.e. are coded in protein
gunzip -c ${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  > $cds_file

# Prep the SNV consensus data for evaluation downstream
Rscript --vanilla ${analysis_dir}/00-tp53-nf1-alterations.R \
  --snvConsensus ${consensus_file} \
  --clinicalFile ${clinical_file} \
  --outputFolder ${analysis_dir}/results \
  --gencode ${cds_file}

# expression files for prediction
collapsed_stranded="pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
collapsed_polya="pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"

# Run classifier for stranded and polya
python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_stranded}
python3 ${analysis_dir}/01-apply-classifier.py -f ${collapsed_polya}


# Run ROC plot step

# Skip poly-A plotting in CI
if [ "$POLYA_PLOT" -gt "0" ]; then
	python3 ${analysis_dir}/02-evaluate-classifier.py -s ${analysis_dir}/results/TP53_NF1_snv_alteration.tsv -f ${analysis_dir}/results/pbta-gene-expression-rsem-fpkm-collapsed.polya_classifier_scores.tsv -c ${data_dir}/pbta-histologies.tsv -o polya
fi

python3 ${analysis_dir}/02-evaluate-classifier.py -s ${analysis_dir}/results/TP53_NF1_snv_alteration.tsv -f ${analysis_dir}/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded_classifier_scores.tsv -c ${data_dir}/pbta-histologies.tsv -o stranded
