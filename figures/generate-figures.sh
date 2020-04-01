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
data_dir="$BASEDIR/data"

# Make output folders for all figures
mkdir -p pngs

################ Sample distribution
# Run sample distribution analysis
bash ${analyses_dir}/sample-distribution-analysis/run-sample-distribution.sh

# Run the figure assembly
Rscript scripts/fig1-sample-distribution.R

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

####### Transcriptomic overview

# First run the dimension reduction steps 
bash ${analyses_dir}/transcriptomic-dimension-reduction/dimension-reduction-plots.sh

# Then collapse RNA-seq data, which is required for GSVA and immune deconvolution
bash ${analyses_dir}/collapse-rnaseq/run-collapse-rnaseq.sh

# Generate GSVA scores
Rscript --vanilla ${analyses_dir}/gene-set-enrichment-analysis/01-conduct-gsea-analysis.R \
  --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
  --output ${analyses_dir}/gene-set-enrichment-analysis/results/gsva_scores_stranded.tsv

# Immune deconvolution - we can't use CIBERSORT because we don't have access to it
Rscript --vanilla ${analyses_dir}/immune-deconv/01-immune-deconv.R \
  --polyaexprs ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds \
  --strandedexprs ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
  --clin ${data_dir}/pbta-histologies.tsv \
  --method "mcp_counter" \
  --output ${analyses_dir}/immune-deconv/results/deconv-output-for-figures.RData
