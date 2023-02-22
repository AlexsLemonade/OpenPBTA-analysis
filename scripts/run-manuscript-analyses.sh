#!/bin/bash
#
# S. Spielman for CCDL ALSF, 2022-2023
# Run all analysis modules that are used in the manuscript

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

# Analyses are run in alphabetical order, except when their results are required
#  as input to other modules. For analysis details, see:
#  https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/README.md

# First, run modules that cannot be run locally due to memory requirements
if [ "$RUN_LOCAL" -lt "1" ]; then

  # Run the `focal-cn-file-preparation` module
  bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh

  # Run the `snv-callers` module
  #  Set OPENPBTA_MANUSCRIPT=1 to overwrite any existing tables, databases
  OPENPBTA_MANUSCRIPT=1 bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-pbta.sh
  OPENPBTA_MANUSCRIPT=1 bash ${analyses_dir}/snv-callers/run_caller_consensus_analysis-tcga.sh
fi

# Next, run modules that can be run locally

# Run the `chromosomal-instability` module
#  Generates input required for `chromothripsis` module
bash ${analyses_dir}/chromosomal-instability/run_breakpoint_analysis.sh

# Run the `chromothripsis` module
bash ${analyses_dir}/chromothripsis/run-chromothripsis.sh

# Run the `cnv-chrom-plot` module, which contains an Rmd file that generates
#  `cn_status_bp_per_bin.tsv` results for use in Figure S3 heatmap
Rscript -e "rmarkdown::render('${analyses_dir}/cnv-chrom-plot/cn_status_heatmap.Rmd')"

# Run the `fusion_filtering` module
bash ${analyses_dir}/fusion_filtering/run_fusion_merged.sh

# Run the `gene-set-enrichment-analysis` module
bash ${analyses_dir}/gene-set-enrichment-analysis/run-gsea.sh

# Run the `immune-deconv` module
#  Generates input required for the `survival-analysis` module
bash ${analyses_dir}/immune-deconv/run-immune-deconv.sh

# Run the `interaction-plots` module
bash ${analyses_dir}/interaction-plots/01-create-interaction-plots.sh

# Run the `tp53_nf1_score` module
#  Generates input required for `mutational-signatures` and `survival-analysis` modules
bash ${analyses_dir}/tp53_nf1_score/run_classifier.sh

# Run the `tp53_nf1_score` module at tumor purity threshold
#  Generates panels and results used to make additional panels for Figure S7
bash ${analyses_dir}/tp53_nf1_score/run_classifier-tumor_purity_threshold.sh

# Run the `mutational-signatures` module
#  This only runs the part of the module used in the manuscript (i.e., not de novo)
OPENPBTA_CNS_FIT_ONLY=1 bash ${analyses_dir}/mutational-signatures/run_mutational_signatures.sh

# Run the `oncoprint-landscape` module
bash ${analyses_dir}/oncoprint-landscape/run-oncoprint.sh

# Run the `run-gistic` module
bash ${analyses_dir}/run-gistic/run-gistic-module.sh

# Run the `survival-analysis` module
bash ${analyses_dir}/survival-analysis/run_survival.sh

# Run the `telomerase-activity-prediction` module
OPENPBTA_FOR_FIGURES=1 bash ${analyses_dir}/telomerase-activity-prediction/RUN-telomerase-activity-prediction.sh

# Run the `transcriptomic-dimension-reduction` module
bash ${analyses_dir}/transcriptomic-dimension-reduction/dimension-reduction-plots.sh

# Run the `tumor-purity-exploration` module
bash ${analyses_dir}/tumor-purity-exploration/run_tumor-purity.sh
