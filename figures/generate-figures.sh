#!/bin/bash
#
# Run all figure making scripts.

#enviroment settings
set -e
set -o pipefail

# Find current directory based on this script
WORKDIR=$(dirname "${BASH_SOURCE[0]}")
cd "$WORKDIR"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"
data_dir="$BASEDIR/data"
scratch_dir="$BASEDIR/scratch"

# Make output folders for all figures
mkdir -p pngs

################ Sample distribution
# Run sample distribution analysis
bash ${analyses_dir}/sample-distribution-analysis/run-sample-distribution.sh

# Run the figure assembly
Rscript --vanilla scripts/fig1-sample-distribution.R

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
Rscript --vanilla -e "rmarkdown::render('../analyses/mutational-signatures/mutational_signatures.Rmd', clean = TRUE)"

# Run the figure assembly
Rscript --vanilla scripts/fig2-mutational-landscape.R

######################
## Interaction plots

# Run the main figure generation script
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh

# Copy the main figure to final directory
cp ${analyses_dir}/interaction-plots/plots/combined_top50.png pngs/mutation_cooccurrence_figure.png

######################
## Oncoprint plot(s)

# Run the `focal-cn-file-preparation` module shell script to prepare the focal
# CN file so that it can be represented on the oncoprint
bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh

# Run the `oncoprint-landscape` module shell script
bash ${analyses_dir}/oncoprint-landscape/run-oncoprint.sh

# Copy the primary plus samples oncoprint figure to final directory
cp ${analyses_dir}/oncoprint-landscape/plots/all_participants_primary-plus_goi_oncoprint.png pngs/fig3-oncoprint-landscape.png

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

# This notebook tests for differences between `short_histology`
# and `integrated_diagnosis` - we can use this to filter what pathways are
# displayed in a heatmap
Rscript --vanilla -e "rmarkdown::render('${analyses_dir}/gene-set-enrichment-analysis/02-model-gsea.Rmd', clean = TRUE)"

# Immune deconvolution - we can't use CIBERSORT because we don't have access to it
# By not supplying an argument to --method, we are electing only to use xCell
Rscript --vanilla ${analyses_dir}/immune-deconv/01-immune-deconv.R \
  --polyaexprs ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds \
  --strandedexprs ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
  --clin ${data_dir}/pbta-histologies.tsv \
  --output ${analyses_dir}/immune-deconv/results/deconv-output-for-figures.RData

# Step that generates the transcriptomic overview figure itself
Rscript --vanilla scripts/transcriptomic-overview.R

####### CN Status Heatmap
# Run consensus CNV so we have a refreshed `pbta-cnv-consensus.seg.gz` file
bash ${analyses_dir}/copy_number_consensus_call/run_consensus_call.sh

# Run CN status heatmap but use parameter so file is saved to figures folder
Rscript -e "rmarkdown::render('${analyses_dir}/cnv-chrom-plot/cn_status_heatmap.Rmd',
                              clean = TRUE, params = list(final_figure=TRUE))"

							  
####### Telomerase Activities


#generate telomerase activities using gene expression data from collapse RNA seq data files
Rscript --vanilla 01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_FPKM.txt
Rscript --vanilla 01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAPolya_FPKM.txt
Rscript --vanilla 01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_counts.txt
Rscript --vanilla 01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAPolya_counts.txt


Rscript --vanilla scripts/TelomeraseActivitites.R

