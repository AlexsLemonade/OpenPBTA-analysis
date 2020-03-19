#!/bin/bash
# 
# Run all figure making scripts. 

# Find current directory based on this script
cd "$(dirname "${BASH_SOURCE[0]}")"
analyses_dir=../analyses

# Make output folders for all figures
mkdir -p pngs
mkdir -p pdfs

## Sample distribution

################ Mutational landscape figure
# Run both SNV caller consensus scripts
# Note: This the PBTA consensus script requires at least 128 MB of RAM to run
bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-pbta.sh
bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-tcga.sh

# Run mutational signatures analysis
Rscript -e "rmarkdown::render('../analyses/mutational-signatures/mutational_signatures.Rmd', clean = TRUE)"

# Run the figure assembly
Rscript scripts/fig2-mutational-landscape.R

######################
## Interaction plots

# Run the main figure generation script
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh

# Copy the main figure to final directory
cp ${analyses_dir}/interaction-plots/plots/combined_top50.pdf pdfs/mutation_cooccurrence_figure.pdf





## Oncoprint plot(s)

## Copy number status heatmap

## Transcriptomic overview
