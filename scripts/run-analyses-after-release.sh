#!/bin/bash
#
# Run all analysis modules that should be re-run after molecular subtyping.
# This script ensures relevant analysis modules are up-to-date following
#  a data release.

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
  # Run the `focal-cn-file-preparation`` module without `OPENPBTA_BASE_SUBTYPING`
  bash ${analyses_dir}/focal-cn-file-preparation/run-prepare-cn.sh
fi

# Run the `fusion_filtering` module without `OPENPBTA_BASE_RELEASE`
bash ${analyses_dir}/fusion_filtering/run_fusion_merged.sh

# Run the `run-gistic` module without `OPENPBTA_BASE_RELEASE`
bash ${analyses_dir}/run-gistic/run-gistic-module.sh
