#!/bin/bash
# 
# Run all figure making scripts. 

# Find current directory based on this script
WORKDIR=$(dirname "${BASH_SOURCE[0]}")
cd "$WORKDIR"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"

# Make output folders for all figures
mkdir -p pngs

## Sample distribution

################ Mutational landscape figure
# Run both SNV caller consensus scripts
# Note: This the PBTA consensus script requires at least 128 MB of RAM to run
# These scripts are intended to run from the base directory, 
# so we will temporarily move there
cd $BASEDIR
bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-pbta.sh
bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-tcga.sh
cd $WORKDIR

# Run mutational signatures analysis
Rscript -e "rmarkdown::render('../analyses/mutational-signatures/mutational_signatures.Rmd', clean = TRUE)"

# Run the figure assembly
Rscript scripts/fig2-mutational-landscape.R

######################
## Interaction plots

# Run the main figure generation script
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh

# Copy the main figure to final directory
cp ${analyses_dir}/interaction-plots/plots/combined_top50.png pngs/mutation_cooccurrence_figure.png





## Oncoprint plot(s)

## Copy number status heatmap

## Transcriptomic overview
