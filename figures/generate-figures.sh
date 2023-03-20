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
mkdir -p pdfs/fig1/panels
mkdir -p pdfs/fig2/panels
mkdir -p pdfs/fig3/panels
mkdir -p pdfs/fig4/panels
mkdir -p pdfs/fig5/panels
mkdir -p pdfs/supp/figs1/panels
mkdir -p pdfs/supp/figs2/panels
mkdir -p pdfs/supp/figs3/panels
mkdir -p pdfs/supp/figs4/panels
mkdir -p pdfs/supp/figs5/panels
mkdir -p pdfs/supp/figs6/panels
mkdir -p pdfs/supp/figs7/panels

#### Make sure color palettes are up-to-date ------------------------
Rscript --vanilla scripts/color_palettes.R
Rscript -e "rmarkdown::render('mapping-histology-labels.Rmd', clean = TRUE, params = list(release = 'release-v23-20230115'))"

### Make sure relevant scratch/ files are up to date ------------------------------

# We will check that necessary input files exist and are <=10 days old,
#  assuming that <=10 days is likely up-to-date.
#   If input files are older, exit with code 0 and prompt to run the module.
MAX_DAYS=10

# Check for an oncoprint file:
oncoprint_file_check=${scratch_dir}/oncoprint_files/primary_only_maf.tsv # use one of the input files in the check
if [ ! -f $oncoprint_file_check ]; then
  echo "The 'oncoprint-landscape' module scratch files do not all exist. Please re-run the 'oncoprint-landscape' module first, or run the script 'scripts/run-manuscript-analyses.sh'."
  exit 1
fi
oncoprint_scratch_days=$((($(date +%s) - $(date +%s -r "$oncoprint_file_check")) / 86400))  # https://unix.stackexchange.com/questions/102691/get-age-of-given-file
if [ $oncoprint_scratch_days -gt ${MAX_DAYS} ]; then
  echo "WARNING: The 'oncoprint-landscape' module scratch files may be out of date. You may want to re-run the 'oncoprint-landscape' module first, or run the script 'scripts/run-manuscript-analyses.sh'."
fi

# Check for snv-callers databases:
snv_pbta_db_check=${scratch_dir}/snv-callers/snv_db.sqlite
snv_tcga_db_check=${scratch_dir}/snv-callers/tcga_snv_db.sqlite
if [ ! -f $snv_pbta_db_check ] || [ ! -f $snv_tcga_db_check ]; then
  echo "The 'snv-callers' module scratch databases do not exist. Please re-run the 'snv-callers' module first, or run the script 'scripts/run-manuscript-analyses.sh'."
  exit 1
fi

snv_pbta_scratch_days=$((($(date +%s) - $(date +%s -r "$snv_pbta_db_check")) / 86400))  # https://unix.stackexchange.com/questions/102691/get-age-of-given-file
snv_tcga_scratch_days=$((($(date +%s) - $(date +%s -r "$snv_tcga_db_check")) / 86400))  # https://unix.stackexchange.com/questions/102691/get-age-of-given-file
if [ $snv_pbta_scratch_days -gt ${MAX_DAYS} ] || [ $snv_tcga_scratch_days -gt ${MAX_DAYS} ]; then
  echo "WARNING: The 'snv-callers' module scratch databases may be out of date. You may want to re-run the 'snv-callers' module first, or run the script 'scripts/run-manuscript-analyses.sh'."
fi


############################# Figure panels ##########################################

# Below, we establish figure panels in `figures/pdf/` for all manuscript figures.
# For each, we first run any `figures/scripts/` scripts that generate figure panels,
#  and then we copy any additional panels that were generated in `analyses/` modules.

##### Figure 1: Workflow and sample distribution ------------------------------------

# Generate sample distribution panel for Figure 1
Rscript --vanilla scripts/fig1-sample-distribution.R


##### Figure 2: Oncoprint ------------------------------------------------------------

# Create single panel PDFs and legends
# Note: This script also generates a panel Figure S3A
Rscript --vanilla scripts/fig2_figS3-oncoprint-landscape.R


##### Figure 3: Mutation overview ----------------------------------------------------

# Copy the main figure to final directory - panels A,B
cp ${analyses_dir}/interaction-plots/plots/combined_top50.pdf pdfs/fig3/panels/mutation_cooccurrence_figure.pdf

# Copy scatter plot from the chromothripsis module - Panel C
cp ${analyses_dir}/chromothripsis/plots/04-breakpoint-data/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf pdfs/fig3/panels/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf

# Run the Rscript that creates the barplot using the most recent color palette - panel D
Rscript --vanilla scripts/fig3-chromothripsis-barplot.R

# Run the Rscript that creates the sina and IQR plot - panel E
# This script also creates both panels for S4
Rscript --vanilla scripts/fig3_figS4-mutational-signatures-panels.R


##### Figure 4: TP53 and telomerase ---------------------------------------------------


# Generate TP53 and telomerase figures
Rscript --vanilla scripts/fig4-tp53-telomerase-panels.R # panels A, B, C, D, F
Rscript --vanilla scripts/fig4-heatmap.R  # panel E
Rscript --vanilla scripts/fig4-hgg-subtype-forest-plot.R # panel G
Rscript --vanilla scripts/fig4-hgg-kaplan-meier.R # panel H


##### Figure 5: GSVA and immune deconvolution -----------------------------------------

# Generate the GSVA, UMAP, and legend panels: panels A and B
Rscript --vanilla scripts/fig5-panels-gsva-umap.R

# Generate panels C and E:
# This script additionally generates panels S6E and S6F
Rscript --vanilla scripts/fig5_figS6-immune-deconv-panels.R

# Generate the forest plot for panel D
Rscript --vanilla scripts/fig5-forest-plot.R




##### Figure S2: Consensus SNV calls and TMB --------------------------------------------

# Generate SNV figures, **but only if NOT LOCAL**. Significant memory requirements.
if [ "$RUN_LOCAL" -lt "1" ]; then
  Rscript --vanilla scripts/figS2-snv-callers-panels.R  # Figure S2 panels A-G
fi
# Generate TMB panels, which does not have the same significant memory requirements
Rscript --vanilla scripts/figS2-tmb-compare-panels.R  # Figure S2 panels H, I






##### Figure S3: Other oncoprint and CNV landscape --------------------------------------

# Copy tumor purity violin across cancer groups plot tumor purity exploration module - Panel A
cp ${analyses_dir}/tumor-purity-exploration/plots/cancer_group_tumor_fraction.pdf pdfs/supp/figs3/panels/cancer_group_tumor_fraction.pdf

# The Panel B oncoprint panel was already generated by a Figure 2 step: `Rscript --vanilla scripts/fig2_figS3-oncoprint-landscape.R`

# Generate remaining three panels for S3 (C, D, E): CN status heatmap and chromothripsis region boxplots
Rscript --vanilla scripts/figS3-panels-cn-chromothripsis.R




##### Figure S4: More mutational signatures --------------------------------------------

# There is no code in this section because Figure S4 panels were previously generated by
#  the script `fig3_figS4-mutational-signatures-panels.R` run during Figure 3 steps



##### Figure S5: More TP53/telomerase --------------------------------------------------


# Generate all figure panels A,B,C
Rscript --vanilla scripts/figS5-all-panels.R




##### Figure S6: More UMAP and other molecular subtypes --------------------------------


# Generate UMAP panels - A, B, C, D
Rscript --vanilla scripts/figS6-subtype-umap-panels.R


# Panels E and F were previously generated with `scripts/fig5_figS6-immune-deconv-panels.R` in a Figure 5 step


##### Figure S7: Transcriptomics results for only high tumor-purity samples -------------

# Run script to generate Panel A: Barplot of cancer groups by RNA library preparation
Rscript --vanilla scripts/figS7-seqcenter-barplot.R

# Run script to generate Panel B: UMAP showing both polyA and stranded RNA library preparations
Rscript --vanilla scripts/figS7-UMAP-libraries.R

# Run script to generate Panel G: TP53 and EXTEND score boxplots and associated legend
Rscript --vanilla scripts/figS7-tp53-telomerase-tumor-purity-threshold.R

# Copy panels for this plot:

# Panel C: UMAP highlight sequencing centers
cp ${analyses_dir}/transcriptomic-dimension-reduction/plots/umap_cancer-group_sequencing-center.pdf pdfs/supp/figs7/panels/


# Panel D: TP53 classifier ROC curve
cp ${analyses_dir}/tp53_nf1_score/results/tumor-purity-threshold/tp53-roc_tumor-purity-threshold.pdf pdfs/supp/figs7/panels/

# Panel E: TP53 scores across TP53 status
cp ${analyses_dir}/tp53_nf1_score/results/tumor-purity-threshold/tp53-score-status_tumor-purity-threshold.pdf pdfs/supp/figs7/panels/

# Panel F: TP53 expression across TP53 status
cp ${analyses_dir}/tp53_nf1_score/results/tumor-purity-threshold/tp53-fpkm-status_tumor-purity-threshold.pdf pdfs/supp/figs7/panels/

# Panel H: UMAP highlighting broad histologies
cp ${analyses_dir}/transcriptomic-dimension-reduction/plots/umap_tumor-purity-threshold.pdf pdfs/supp/figs7/panels/

# Panel I: quanTIseq fractions across cancer groups
cp ${analyses_dir}/immune-deconv/plots/tumor-purity-threshold_quantiseq-cancer-groups.pdf pdfs/supp/figs7/panels/



##### Clean up --------------------------------------------------------------------------

# Sometimes Rplots.pdf gets produced and honestly nobody really knows why.
rm -f Rplots.pdf # use `-f` to not get a warning if Rplots.pdf is NOT there.

