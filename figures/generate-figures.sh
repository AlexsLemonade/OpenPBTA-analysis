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


##### Figure 1: Workflow and sample distribution ---------------------------

# Create directories
mkdir -p pdfs/fig1/panels
mkdir -p pdfs/figs1/panels

# Generate sample distribution panel for Figure 1 (and supplementary panels)
Rscript --vanilla scripts/fig1-sample-distribution.R


##### Figure 2: Oncoprint ---------------------------

# Create directory
mkdir -p pdfs/fig2/panels

if [ "$RUN_LOCAL" -lt "1" ]; then
  # Run the `focal-cn-file-preparation` module shell script to prepare the focal
  # CN file so that it can be represented on the oncoprint
  bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh
fi

# Run the `oncoprint-landscape` module shell script
bash ${analyses_dir}/oncoprint-landscape/run-oncoprint.sh


# Create single panel PDFs and legends
# Note: This script also generates a panel for Figure S3
Rscript --vanilla scripts/fig2-oncoprint-landscape.R

##### Figure 3: Mutation overview ---------------------------

#### Interaction plots #####

# Create directory
mkdir -p pdfs/fig3/panels

# Run the main figure generation script
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh

# Copy the main figure to final directory
cp ${analyses_dir}/interaction-plots/plots/combined_top50.pdf pdfs/fig3/panels/mutation_cooccurrence_figure.pdf


#### Chromothripsis ####
# The chromothripsis module uses breakpoint counts from this module
bash ${analyses_dir}/chromosomal-instability/run_breakpoint_analysis.sh

# Run the chromothripsis module
bash ${analyses_dir}/chromothripsis/run-chromothripsis.sh

# Copy scatter plot from the chromothripsis module
cp ${analyses_dir}/chromothripsis/plots/04-breakpoint-data/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf pdfs/fig3/panels/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf

# Run the Rscript that creates the barplot using the most recent color palette
Rscript --vanilla scripts/fig3-chromothripsis-barplot.R

#### Mutational signatures ####

# Run the module
bash ${analyses_dir}/mutational/run_mutational_signatures.sh


# Copy panel from module
cp ${analyses_dir}/mutational-signatures/plots/cns/exposures_sina_IQR.pdf  pdfs/fig3/panels/mutational_signatures_exposures.pdf


##### Figure 4: TP53 and telomerase ---------------------------

# Create directory
mkdir -p pdfs/fig4/panels

#### Generate files for telomerase figures ####

# Generate collapsed data for count files
Rscript ${analyses_dir}/collapse-rnaseq/01-summarize_matrices.R \
  -i ${data_dir}/pbta-gene-counts-rsem-expected_count.stranded.rds \
  -g ${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz \
  -m ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds \
  -t ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.stranded.rds

# Generate telomerase activities using gene expression data from collapse RNA seq data files
Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_FPKM.txt

#### Generate files for tp53 figures ####
bash ${analyses_dir}/tp53_nf1_score/run_classifier.sh


# Generate TP53 and telomerase figures 
Rscript --vanilla scripts/fig4-tp53-telomerase.R



##### Figure 5: GSVA and immune deconvolution ---------------------------

# Create directory
mkdir -p pdfs/fig5/panels


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


##### Figure S2: Consensus SNV calls and TMB ---------------------------
##### CAUTION!!! Generating these panels requires 32 GB RAM!

# Run the snv module for PBTA and TCGA
bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-pbta.sh
bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-tcga.sh

# Generate SNV and TMB figures
Rscript --vanilla scripts/supp-snv-callers-panels.R
Rscript --vanilla scripts/supp-tmb-compare-panels.R



##### Figure S3: Other oncoprint and CNV landscape ---------------------------


# Create directory
mkdir -p pdfs/supp/figs3/panels


# The oncoprint panel for this figure was already generated by a Figure 2 step: `Rscript --vanilla scripts/fig2-oncoprint-landscape.R`


#### Prepare for CN status heatmap figure ####
if [ "$RUN_LOCAL" -lt "1" ]; then
    # Run consensus CNV so we have a refreshed `pbta-cnv-consensus.seg.gz` file
    bash ${analyses_dir}/copy_number_consensus_call/run_consensus_call.sh
fi

# Run CN status heatmap WITHOUT final file parameter, to generate `cn_status_bp_per_bin.tsv` input file for figure scripts
Rscript -e "rmarkdown::render('${analyses_dir}/cnv-chrom-plot/cn_status_heatmap.Rmd', clean = TRUE)"

# Chromothripsis data was already generated by a Figure 3 step: `bash ${analyses_dir}/chromosomal-instability/run_breakpoint_analysis.sh`

 
# Generate three panels for S3
Rscript --vanilla scripts/supp-S3-panels-BCD.R

##### Figure S4: More mutational signatures ---------------------------

# Create directory
mkdir -p pdfs/supp/figs4/panels

# The mutational signatures data (and plots) were already generated by a Figure 3 step: `bash ${analyses_dir}/mutational-signatures/run_mutational_signatures.sh`

# Copy panels to directory
cp ${analyses_dir}/mutational-signatures/plots/cns/exposures_per_sample_barplot.pdf  pdfs/supp/figs4/panels/
cp ${analyses_dir}/mutational-signatures/plots/cns/signature1_tumor-descriptor_cancer-groups.pdf  pdfs/supp/figs4/panels/



##### Figure S5: Other TP53, telomerase, and UMAPs ---------------------------

# Create directory
mkdir -p pdfs/supp/figs5/panels

# The TP53 data needed for ROC was already generated by a Figure 4 step: `bash ${analyses_dir}/tp53_nf1_score/run_classifier.sh`
# The FPKM telomerase data needed for scatterplots was already generated by a Figure 4 step: `Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_FPKM.txt`

# Generate first three figure panels
Rscript --vanilla scripts/supp-S5-panels-ABC.R

# The data needed for UMAP plots was already generated by a Figure 5 step: `bash ${analyses_dir}/transcriptomic-dimension-reduction/dimension-reduction-plots.sh`

# Generate UMAP panels
Rscript --vanilla scripts/supp-subtype-umap.R


### Additional figure scripts not yet stubbed out to be in figures:
# scripts/supp-tp53-correlation.R

### Additional preparation code that is not currently used in figures:
###### polya telomerase ########
#Rscript ${analyses_dir}/collapse-rnaseq/01-summarize_matrices.R \
#  -i ${data_dir}/pbta-gene-counts-rsem-expected_count.polya.rds \
#  -g ${data_dir}/gencode.v27.primary_assembly.annotation.gtf.gz \
#  -m ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds \
#  -t ${analyses_dir}/collapse-rnaseq/pbta-gene-counts-rsem-expected_count-collapsed_table.polya.rds
# Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAPolya_FPKM.txt
# Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAPolya_counts.txt




