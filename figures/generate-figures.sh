#!/bin/bash
#
# Run all figure making scripts.

#enviroment settings
set -e
set -o pipefail

# If RUN_LOCAL is used, the snv callers steps are skipped because they cannot
# be run on a local computer -- the idea is that setting RUN_LOCAL=1 will allow for
# local testing running/testing of all the other figures
RUN_LOCAL=${RUN_LOCAL:-0}

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
mkdir -p pdfs

#### Make sure color palettes are up-to-date
Rscript --vanilla scripts/color_palettes.R
Rscript -e "rmarkdown::render('mapping-histology-labels.Rmd', clean = TRUE, params = list(release = 'release-v21-20210820'))"

################ Sample distribution
# Run sample distribution analysis
bash ${analyses_dir}/sample-distribution-analysis/run-sample-distribution.sh

# TODO: Rscript for figure (related: https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/1175)

######################
## Interaction plots

# Create directory that will hold the relevant panel
mkdir -p pdfs/fig3/panels

# Run the main figure generation script
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh

# Copy the main figure to final directory
cp ${analyses_dir}/interaction-plots/plots/combined_top50.pdf pdfs/fig3/panels/mutation_cooccurrence_figure.pdf


#####################
## Chromothripsis

# The chromothripsis module uses breakpoint counts from this module
bash ${analyses_dir}/chromosomal-instability/run_breakpoint_analysis.sh

# Run the chromothripsis module
bash ${analyses_dir}/chromothripsis/run-chromothripsis.sh

# Create directory that will hold the relevant scatter plot from the chromothripsis module
mkdir -p pdfs/fig3/panels
cp ${analyses_dir}/chromothripsis/plots/04-breakpoint-data/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf pdfs/fig3/panels/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf

# Run the Rscript that creates the barplot using the most recent color palette
Rscript --vanilla scripts/fig3-chromothripsis-barplot.R

######################
## Oncoprint plot(s)

if [ "$RUN_LOCAL" -lt "1" ]; then
  # Run the `focal-cn-file-preparation` module shell script to prepare the focal
  # CN file so that it can be represented on the oncoprint
  bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh
fi

# Run the `oncoprint-landscape` module shell script
bash ${analyses_dir}/oncoprint-landscape/run-oncoprint.sh

# Will create two plots - primary only and "primary plus" samples
filenames=(primary_only primary-plus)

# Create single panel PDFs and legends
Rscript --vanilla scripts/fig2-oncoprint-landscape.R

####### Telomerase Activities

# Generate collapsed data for count files
Rscript ${analyses_dir}/collapse-rnaseq/01-summarize_matrices.R \
  -i ${data_dir}/pbta-gene-counts-rsem-expected_count.stranded.rds \
  -g ${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz \
  -m ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds \
  -t ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.stranded.rds

Rscript ${analyses_dir}/collapse-rnaseq/01-summarize_matrices.R \
  -i ${data_dir}/pbta-gene-counts-rsem-expected_count.polya.rds \
  -g ${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz \
  -m ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds \
  -t ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.polya.rds

# Generate telomerase activities using gene expression data from collapse RNA seq data files
Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_FPKM.txt
Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAPolya_FPKM.txt
Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_counts.txt
Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAPolya_counts.txt

# Build figures of telomerase activity
Rscript --vanilla scripts/fig4-telomerase-activities.R

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

# Step that generates the GSVA, UMAP, and legend panels
Rscript --vanilla scripts/fig5-panels-gsva-umap.R


###### TP53 scores

# Generate the tp53 scores boxplot for figure 3
Rscript --vanilla scripts/fig3-panel-tp53.R

###### Hypermutator signatures
# Run the mutational-signatures module
bash ${analyses_dir}/mutational-signatures/run_mutational_signatures.sh


# Copy the figure to final directory
cp ${analyses_dir}/mutational-signatures/plots/cns/hypermutator_sigs_heatmap.pdf  pdfs/fig4/panels/hypermutator_sigs_heatmap.pdf
cp ${analyses_dir}/mutational-signatures/plots/cns/hypermutator_sigs_heatmap_legends.pdf  pdfs/fig4/panels/hypermutator_sigs_heatmap_legends.pdf

####### Sample distributions

# Generate sample distribution panel for Figure 1 and supplementary panels
Rscript --vanilla scripts/fig1-sample-distribution.R


####### CN Status Heatmap
if [ "$RUN_LOCAL" -lt "1" ]; then
# Run consensus CNV so we have a refreshed `pbta-cnv-consensus.seg.gz` file
bash ${analyses_dir}/copy_number_consensus_call/run_consensus_call.sh
fi

# Run CN status heatmap but use parameter so file is saved to figures folder
Rscript -e "rmarkdown::render('${analyses_dir}/cnv-chrom-plot/cn_status_heatmap.Rmd',
                              clean = TRUE, params = list(final_figure=TRUE))"



######## Mutational signatures
# Copy the figure to final directory
cp ${analyses_dir}/mutational-signatures/plots/cns/exposures_sina_IQR.pdf  pdfs/fig3/panels/mutational_signatures_exposures.pdf


######## Immune deconvolution with quanTIseq (5C)
# run the immune-deconv module:
bash ${analyses_dir}/immune-deconv/run-immune-deconv.sh
# copy figure panel:
cp ${analyses_dir}/immune-deconv/plots/cell_types-cancer_groups.pdf  pdfs/fig5/panels/quantiseq-cell_types-cancer_groups.pdf

###### Survival plots for 4D
Rscript --vanilla scripts/fig4-hgg-kaplan-meier.R
Rscript --vanilla scripts/fig4-hgg-subtype-forest-plot.R
Rscript --vanilla scripts/fig4-tp53-telomerase-panels.R
Rscript --vanilla scripts/fig5-forest-plot.R

###### Box plot for 5E
cp ${analyses_dir}/immune-deconv/plots/cd274_expression_mb_subtypes.pdf  pdfs/fig5/panels/cd274_expression_mb_subtypes.pdf


####### Supplementary figures


# Panels for supplementary figure 3
Rscript --vanilla scripts/supp-S3-panels-BCD.R


# Copy Figure S4 panels (analysis module was run previously)
cp ${analyses_dir}/mutational-signatures/plots/cns/signature1_tumor-descriptor_cancer-groups.pdf   pdfs/supp/figs4/panels/
cp ${analyses_dir}/mutational-signatures/plots/cns/exposures_per_sample_barplot.pdf    pdfs/supp/figs4/panels/


# UMAP panels for supplementary figure 6 from molecular analysis 
Rscript --vanilla scripts/supp-subtype-umap.R

# Copy additional S6 panels (analysis module was run previously)
# 6S - E
cp ${analyses_dir}/immune-deconv/plots/cell_types-molecular_subtypes.pdf pdfs/supp/figs6/panels/quantiseq-cell_types-molecular_subtypes.pdf

# 6S - F
cp ${analyses_dir}/immune-deconv/plots/cd8_cd4_ratio.pdf pdfs/supp/figs6/panels/cd8_cd4_ratio.pdf



