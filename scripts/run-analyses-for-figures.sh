#!/bin/bash
#
# Run all analysis modules that are used for manuscript figures

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

#################################

# Run modules that cannot be run locally due to memory requirements
if [ "$RUN_LOCAL" -lt "1" ]; then
  # This module is not strictly needed for figure scripts, but we should ensure it's up-to-date
  #  and re-run after subtyping. That can be done here elsewhere; commenting out for now.
  #bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh

  # Must be run to SQL databases needed to make Figure S2 panels
  bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-pbta.sh # Figure S2
  bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-tcga.sh # Figure S2
fi

# Run the `oncoprint-landscape` module shell script, for Figure 2 and S3
bash ${analyses_dir}/oncoprint-landscape/run-oncoprint.sh

# Run the interaction plots script for Figures 3 and S3
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh


# Run the chromothripsis module, which requires the chromosomal-instability module, for Figure 3
bash ${analyses_dir}/chromosomal-instability/run_breakpoint_analysis.sh
bash ${analyses_dir}/chromothripsis/run-chromothripsis.sh

# Run the tp53 classifier, for Figures 4 and S5
#  Note this module is run earlier in this script than the relevant
#  figures' placement because these modules use the TP53 scores:
#  `mutational-signatures` and `survival-analysis`
bash ${analyses_dir}/tp53_nf1_score/run_classifier.sh


# Run the mutational-signatures module for Figures 3 and S4
# We only run the part of the module used in the manuscript (i.e., not de novo)
OPENPBTA_CNS_FIT_ONLY=1 bash ${analyses_dir}/mutational-signatures/run_mutational_signatures.sh

# Run the telomerase activity prediction script, for Figures 4 and S5
OPENPBTA_FOR_FIGURES=1 bash ${analyses_dir}/telomerase-activity-prediction/RUN-telomerase-activity-prediction.sh

# Run the immune deconvolution module, for Figures 5 and S6
#  Note this module is run earlier in this script than the relevant
#  figures' placement because the module `survival-analysis`
#  uses the deconvolution results
bash ${analyses_dir}/immune-deconv/run-immune-deconv.sh

# Run the survival module, for Figures 4 and 5
bash ${analyses_dir}/survival-analysis/run_survival.sh

# Run the dimension reduction module, for Figures 5 and S6
bash ${analyses_dir}/transcriptomic-dimension-reduction/dimension-reduction-plots.sh

# Generate GSVA scores and test for cancer group differences, for Figure 5
bash ${analyses_dir}/gene-set-enrichment-analysis/run-gsea.sh

# Generate `cn_status_bp_per_bin.tsv` results for Figure S3 heatmap
Rscript -e "rmarkdown::render('${analyses_dir}/cnv-chrom-plot/cn_status_heatmap.Rmd')"


