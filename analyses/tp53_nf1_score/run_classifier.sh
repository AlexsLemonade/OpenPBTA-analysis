#!/bin/bash

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

data_dir="../../data"
scratch_dir="../../scratch"
# cds gencode bed file  
exon_file="${scratch_dir}/gencode.v27.primary_assembly.annotation.bed"
consensus_file="${data_dir}/pbta-snv-consensus-mutation.maf.tsv.gz"
clinical_file="${data_dir}/pbta-histologies.tsv"
# expression files for prediction
collapsed_stranded="${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"
collapsed_polya="${data_dir}/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"

# Convert GTF to BED file
# Here we are only extracting lines with as a CDS i.e. are coded in protein
gunzip -c ${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz \
  | awk '$3 ~ /CDS/' \
  | convert2bed --do-not-sort --input=gtf - \
  > $exon_file

# Prep the SNV consensus data for evaluation downstream
Rscript --vanilla 00-tp53-nf1-alterations.R \
  --snvConsensus ${consensus_file} \
  --clinicalFile ${clinical_file} \
  --outputFolder results \
  --gencode ${exon_file}

# Run classifier for stranded and polya
python3 01-apply-classifier.py -f ${collapsed_stranded}
python3 01-apply-classifier.py -f ${collapsed_polya}

