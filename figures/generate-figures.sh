#!/bin/bash
#
# Run all figure making scripts.

#enviroment settings
set -e
set -o pipefail

# If RUN_LOCAL is used, the time-intensive steps are skipped because they cannot
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


##### Run all analysis modules required for figures ---------------------------------------

# Run modules that cannot be run locally due to memory requirements
if [ "$RUN_LOCAL" -lt "1" ]; then
  bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh       # Figures 2 and S3
  bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-pbta.sh # Figure S2
  bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-tcga.sh # Figure S2
  bash ${analyses_dir}/copy_number_consensus_call/run_consensus_call.sh  # Figure S3 (heatmap)
fi

# Run the `oncoprint-landscape` module shell script, for Figure 2 and S3
bash ${analyses_dir}/oncoprint-landscape/run-oncoprint.sh

# Run the interaction plots script for Figures 3 and S3
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh


# Run the chromothripsis module, which requires the chromosomal-instability module, for Figure 3
bash ${analyses_dir}/chromosomal-instability/run_breakpoint_analysis.sh
bash ${analyses_dir}/chromothripsis/run-chromothripsis.sh

# Run the mutational-signatures module for Figures 3 and S4
# We only run the part of the module used in the manuscript (i.e., not de novo)
OPENPBTA_CNS_FIT_ONLY=1 bash ${analyses_dir}/mutational-signatures/run_mutational_signatures.sh

# Run the collapse-rnaseq module, which is needed for telomerase, immune deconvolution, and GSVA modules
OPENPBTA_TP53_FIGURES=1 bash ${analyses_dir}/collapse-rnaseq/run-collapse-rnaseq.sh 

# Run the telomerase activity prediction script, for Figures 4 and S5
# TODO: should we actually just run the full module script? I don't do that here since it also re-runs collapse rna seq.
Rscript --vanilla ${analyses_dir}/telomerase-activity-prediction/01-run-EXTEND.R \
 --input ${analyses_dir}/collapse-rnaseq/results/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds \
 --output ${analyses_dir}/telomerase-activity-prediction/results/TelomeraseScores_PTBAStranded_FPKM.txt

# Run the tp53 classifier, for Figures 4 and S5
bash ${analyses_dir}/tp53_nf1_score/run_classifier.sh

# Run the survival module, for Figures 4 and 5
bash ${analyses_dir}/survival-analysis/run_survival.sh

# Run the dimension reduction module, for Figures 5 and S6
bash ${analyses_dir}/transcriptomic-dimension-reduction/dimension-reduction-plots.sh

# Run the immune deconvolution module, for Figures 5 and S6
bash ${analyses_dir}/immune-deconv/run-immune-deconv.sh

# Generate GSVA scores and test for cancer group differences, for Figure 5
bash ${analyses_dir}/gene-set-enrichment-analysis/run-gsea.sh




##### Figure 1: Workflow and sample distribution ------------------------------------

# Create directories
mkdir -p pdfs/fig1/panels

# Generate sample distribution panel for Figure 1 (and supplementary panels)
Rscript --vanilla scripts/fig1-sample-distribution.R


##### Figure 2: Oncoprint ------------------------------------------------------------

# Create directory
mkdir -p pdfs/fig2/panels


# Create single panel PDFs and legends
# Note: This script also generates a panel for Figure S3
Rscript --vanilla scripts/fig2-oncoprint-landscape.R




##### Figure 3: Mutation overview ----------------------------------------------------

# Create directory
mkdir -p pdfs/fig3/panels

# Copy the main figure to final directory - panels A,B
cp ${analyses_dir}/interaction-plots/plots/combined_top50.pdf pdfs/fig3/panels/mutation_cooccurrence_figure.pdf

# Copy scatter plot from the chromothripsis module - Panel C
cp ${analyses_dir}/chromothripsis/plots/04-breakpoint-data/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf pdfs/fig3/panels/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf

# Run the Rscript that creates the barplot using the most recent color palette - panel D
Rscript --vanilla scripts/fig3-chromothripsis-barplot.R

# Copy panel from module - panel E
cp ${analyses_dir}/mutational-signatures/plots/cns/exposures_sina_IQR.pdf  pdfs/fig3/panels/mutational_signatures_exposures.pdf




##### Figure 4: TP53 and telomerase ---------------------------------------------------

# Create directory
mkdir -p pdfs/fig4/panels


# Generate TP53 and telomerase figures 
Rscript --vanilla scripts/fig4-tp53-telomerase-panels.R # panels A, B, C, D, F
Rscript --vanilla scripts/fig4-hgg-subtype-forest-plot.R # panel G
Rscript --vanilla scripts/fig4-hgg-kaplan-meier.R # panel H

# Copy hypermutator figure and legend panel - panel E
cp ${analyses_dir}/mutational-signatures/plots/cns/hypermutator_sigs_heatmap.pdf  pdfs/fig4/panels/hypermutator_sigs_heatmap.pdf
cp ${analyses_dir}/mutational-signatures/plots/cns/hypermutator_sigs_heatmap_legends.pdf  pdfs/fig4/panels/hypermutator_sigs_heatmap_legends.pdf





##### Figure 5: GSVA and immune deconvolution -----------------------------------------

# Create directory
mkdir -p pdfs/fig5/panels


# Generate the GSVA, UMAP, and legend panels: panels A and B
Rscript --vanilla scripts/fig5-panels-gsva-umap.R


# Copy the figures to final directory - panels C and E, in order:
cp ${analyses_dir}/immune-deconv/plots/cell_types-cancer_groups.pdf  pdfs/fig5/panels/quantiseq-cell_types-cancer_groups.pdf
cp ${analyses_dir}/immune-deconv/plots/cd274_expression_mb_subtypes.pdf  pdfs/fig5/panels/cd274_expression_mb_subtypes.pdf


# Generate the forest plot for panel D
Rscript --vanilla scripts/fig5-forest-plot.R







##### Figure S2: Consensus SNV calls and TMB --------------------------------------------

# Create directory
mkdir -p pdfs/supp/figs2/panels

# Generate SNV and TMB figures, **but only if NOT LOCAL**. Signficant memory requirements.
if [ "$RUN_LOCAL" -lt "1" ]; then
  Rscript --vanilla scripts/supp-snv-callers-panels.R   # Figure S2
  Rscript --vanilla scripts/supp-tmb-compare-panels.R   # Figure S2
fi





##### Figure S3: Other oncoprint and CNV landscape --------------------------------------


# Create directory
mkdir -p pdfs/supp/figs3/panels


# Note the oncoprint panel A for this figure was already generated by a Figure 2 step: `Rscript --vanilla scripts/fig2-oncoprint-landscape.R`

# Generate remaining three panels for S3: B, C, and D
Rscript --vanilla scripts/supp-S3-panels-BCD.R






##### Figure S4: More mutational signatures --------------------------------------------

# Create directory
mkdir -p pdfs/supp/figs4/panels

# Copy panel A:
cp ${analyses_dir}/mutational-signatures/plots/cns/exposures_per_sample_barplot.pdf  pdfs/supp/figs4/panels/
# Copy panel B:
cp ${analyses_dir}/mutational-signatures/plots/cns/signature1_tumor-descriptor_cancer-groups.pdf  pdfs/supp/figs4/panels/





##### Figure S5: More TP53/telomerase --------------------------------------------------


# Create directory
mkdir -p pdfs/supp/figs5/panels

# Generate figure panels 
Rscript --vanilla scripts/supp-S5-panels.R






##### Figure S6: More UMAP and other molecular subtypes --------------------------------

# Create directory
mkdir -p pdfs/supp/figs6/panels


# Generate UMAP panels - A, B, C, D
Rscript --vanilla scripts/supp-subtype-umap.R


# Copy panel E
cp ${analyses_dir}/immune-deconv/plots/cell_types-molecular_subtypes.pdf pdfs/supp/figs6/panels/quantiseq-cell_types-molecular_subtypes.pdf
# Copy panel F
cp ${analyses_dir}/immune-deconv/plots/cd8_cd4_ratio.pdf pdfs/supp/figs6/panels/cd8_cd4_ratio.pdf






##### Clean up --------------------------------------------------------------------------

# Sometimes Rplots.pdf gets produced and honestly nobody really knows why.
rm -f Rplots.pdf # use `-f` to not get a warning if Rplots.pdf is NOT there.

