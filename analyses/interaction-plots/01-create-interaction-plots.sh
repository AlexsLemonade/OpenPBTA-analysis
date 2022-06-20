#!/bin/bash

# JA Shapiro for CCDL 2019-2020
#
# Runs scripts/01-process_mutations.R with some default settings.
# Takes one enviroment variable, `OPENPBTA_ALL`, which if 0 runs only
# the full dataset and the largest disease set (for testing). If 1 or more,
# all samples ar run (this is also the default behavior if unset)

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

base_dir=../..
script_dir=scripts
results_dir=results
plot_dir=plots
temp_dir=${base_dir}/scratch/interaction

mkdir -p $temp_dir

ALL=${OPENPBTA_ALL:-1}

ind_samples=${base_dir}/data/independent-specimens.wgs.primary-plus.tsv
metadata=${base_dir}/data/pbta-histologies.tsv

# using consensus
maf=${base_dir}/data/pbta-snv-consensus-mutation.maf.tsv.gz

cooccur=${results_dir}/cooccur_top50
gene_disease=${results_dir}/gene_disease_top50
plot=${plot_dir}/cooccur_top50
disease_plot=${plot_dir}/gene_disease_top50
combined_plot=${plot_dir}/combined_top50

# associative array of diseases to test; chosen by those that are most common
# in the openPBTA dataset
declare -A disease
disease[Diffuse-astrocytic-and-oligodendroglial-tumor]="Diffuse astrocytic and oligodendroglial tumor"
if [ "$ALL" -gt "0" ]; then
  disease[LGAT]="Low-grade astrocytic tumor"
  disease[Tumors-of-sellar-region]="Tumors of sellar region"
  disease[Embryonal-tumor]="Embryonal tumor"
  disease[Ependymal-tumor]="Ependymal tumor"
  disease[Meningioma]="Meningioma"
  disease[Neuronal-mixed-neuronal-glial-tumor]="Neuronal and mixed neuronal-glial tumor"
fi


# Get FLAG file and add header
# include top 50 frequently mutated
exclude_file=FLAGS.tsv
echo gene$'\t'count > $exclude_file
head -n 50 <(curl -s https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5706417/bin/12920_2017_309_MOESM3_ESM.txt)\
  >> $exclude_file


# make output directories if they don't exist
mkdir -p $results_dir
mkdir -p $plot_dir

# run scripts

# all samples first
echo "All"
Rscript ${script_dir}/01-disease-specimen-lists.R \
  --metadata ${metadata} \
  --specimen_list ${ind_samples} \
  --disease "All" \
  --outfile ${temp_dir}/ALL.tsv

Rscript ${script_dir}/02-process_mutations.R \
  --maf ${maf} \
  --metadata ${metadata} \
  --specimen_list ${temp_dir}/ALL.tsv \
  --exclude_genes $exclude_file \
  --vaf 0.05 \
  --min_mutated 5 \
  --max_genes 50 \
  --out ${cooccur}.ALL.tsv \
  --disease_table ${gene_disease}.tsv

Rscript ${script_dir}/03-plot_interactions.R \
  --infile ${cooccur}.ALL.tsv \
  --outfile ${plot}.ALL.pdf \
  --disease_table ${gene_disease}.tsv \
  --disease_plot ${disease_plot}.pdf \
  --combined_plot ${combined_plot}.pdf \
  --plotsize 50

# now individual diseases
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
    --vaf 0.05 \
    --min_mutated 5 \
    --max_genes 50 \
    --out ${cooccur}.${disease_id}.tsv

  Rscript ${script_dir}/03-plot_interactions.R \
    --infile ${cooccur}.${disease_id}.tsv \
    --outfile ${plot}.${disease_id}.png \
    --plotsize 50
done
