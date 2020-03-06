#!/bin/bash
# 
# Run all figure making scripts. 

# Find current directory based on this script
cd "$(dirname "${BASH_SOURCE[0]}")" 
 
# Make an output folder for all directories
mkdir figures

## Sample distribution

################ Mutational landscape figure
# Run both SNV caller consensus scripts
# Note: This the PBTA consensus script requires at least 128 MB of RAM to run
bash analyses/snv-callers/run_caller_consensus_analysis-pbta.sh
bash analyses/snv-callers/run_caller_consensus_analysis-tcga.sh

# Run mutational signatures analysis
Rscript -e "rmarkdown::render('analyses/mutational-signatures/mutational_signatures.Rmd', clean = TRUE)"

# Run the figure assembly
Rscript figure-scripts/fig2-mutational-landscape.R

## Interaction plots 

## Oncoprint plot(s)

## Copy number status heatmap

## Transcriptomic overview
