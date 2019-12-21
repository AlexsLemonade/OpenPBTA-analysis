#!/bin/bash

# JA Shapiro for CCDL 2019
#
# Runs scripts/01-process_mutations.R with some default settings.
# Takes one enviroment variable, `OPENPBTA_ALL`, which if 0 runs only
# the full dataset and the largest disease set (for testing). If 1 or more,
# all samples ar run (this is also the default behavior if unset)

set -e
set -o pipefail

base_dir=analyses/interaction-plots
script_dir=${base_dir}/scripts
results_dir=${base_dir}/results
plot_dir=${base_dir}/plots
temp_dir=scratch

ALL=${OPENPBTA_ALL:-1}

ind_samples=data/independent-specimens.wgs.primary-plus.tsv
metadata=data/pbta-histologies.tsv

# using consensus
maf=data/pbta-snv-consensus-mutation.maf.tsv.gz

cooccur=${results_dir}/consensus_top50
plot=${plot_dir}/consensus_top50

# associative array of diseases to test; chosen by those that are most common
# in the openPBTA dataset
declare -A disease
disease[All]="All"
disease[LGAT]="Low-grade astrocytic tumor"
if [ "$ALL" -gt "0" ]; then
  disease[Medulloblastoma]="Medulloblastoma"
  disease[Ependymoma]="Ependymoma"
  disease[HGAT]="High-grade glioma"
  disease[DIPG]="Brainstem glioma- Diffuse intrinsic pontine glioma"
  disease[Ganglioglioma]="Ganglioglioma"
  disease[Craniopharyngioma]="Craniopharyngioma"
  disease[ATRT]="Atypical Teratoid Rhabdoid Tumor"
fi


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
    --maf ${maf} \
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
