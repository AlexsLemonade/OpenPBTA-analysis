#!/bin/bash

# JA Shapiro for CCDL 2019
#
# Runs scripts/01-process_mutations.R with some default settings.

set -e
set -o pipefail

base_dir=analyses/interaction-plots
script_dir=${base_dir}/scripts
results_dir=${base_dir}/results
plot_dir=${base_dir}/plots
temp_dir=scratch


ind_samples=data/independent-specimens.wgs.primary-plus.tsv
metadata=data/pbta-histologies.tsv

# using lancet data for now
maf=data/pbta-snv-lancet.vep.maf.gz

cooccur=${results_dir}/lancet_top50
plot=${plot_dir}/lancet_top50

declare -A disease
disease[All]="All"
disease[LGAT]="Low-grade glioma;astrocytoma (WHO grade I/II)"
disease[Medulloblastoma]="Medulloblastoma"
disease[Ependymoma]="Ependymoma"
disease[HGAT]="High-grade glioma;astrocytoma (WHO grade III/IV)"
disease[DIPG]="Brainstem glioma- Diffuse intrinsic pontine glioma"
disease[Ganglioglioma]="Ganglioglioma"
disease[Craniopharyngioma]="Craniopharyngioma"
disease[ATRT]="Atypical Teratoid Rhabdoid Tumor"


# make output directories if they don't exist
mkdir -p $results_dir
mkdir -p $plot_dir

# run scripts


for disease_id in "${!disease[@]}"; do
  echo $disease_id
  Rscript ${script_dir}/01-disease-specimen-lists.R \
    --metadata ${metadata} \
    --specimen_list ${ind_samples} \
    --disease "${disease[$disease_id]}" \
    --outfile ${temp_dir}/${disease_id}.tsv

  Rscript ${script_dir}/02-process_mutations.R \
    --maf data/pbta-snv-lancet.vep.maf.gz \
    --metadata ${metadata} \
    --specimen_list ${temp_dir}/${disease_id}.tsv \
    --vaf 0.2 \
    --min_mutated 5 \
    --max_genes 50 \
    --out ${cooccur}.${disease_id}.tsv
  
  Rscript ${script_dir}/03-plot_interactions.R \
    --infile ${cooccur}.${disease_id}.tsv \
    --outfile ${plot}.${disease_id}.png \
    --plotsize 50
done