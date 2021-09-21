#!/bin/bash

# Jaclyn Taroni for ALSF CCDL 2020

# This script should always run as if it were being called from
# the directory it lives in.
cd "$(dirname "${BASH_SOURCE[0]}")"

# If this is CI, run the example included with GISTIC
# The sample size for the subset files are too small otherwise
IS_CI=${OPENPBTA_CI:-0}

if [[ "$IS_CI" -gt "0" ]]
then

# Environmental variables for MCR
ORIG_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/mcr/v83/runtime/glnxa64:/opt/mcr/v83/bin/glnxa64:/opt/mcr/v83/sys/os/glnxa64
export XAPPLRESDIR=/opt/mcr/v83/X11/app-defaults

# We want this to fail if the GISTIC example fails only -- because we have
# some instances of running GISTIC that do not complete but do save some 
# output
set -e
set -o pipefail
# Run the example that comes with GISTIC - that allows us to 
cd /home/rstudio/gistic_install && ./run_gistic_example

# 'Undo' environmental variables for MCR
export LD_LIBRARY_PATH=$ORIG_LD_LIBRARY_PATH
unset XAPPLRESDIR

else

  RESULTSDIR=results
  DATADIR=../../data
  
  # Generate files that are compatible for GISTIC 
  Rscript scripts/prepare_seg_for_gistic.R \
    --in_consensus $DATADIR/cnv-consensus.seg.gz \
    --out_consensus $RESULTSDIR/cnv-consensus-gistic-only.seg.gz \
    --histology $DATADIR/histologies.tsv
  
  # run GISTIC for the whole cohort
  echo "Running GISTIC on the entire OpenPedCan cohort..."
  bash scripts/run-gistic-opentargets.sh

fi

# Since the short_histology field in OT no longer code disease as `LGAT` or `HGAT`, we cannot use that as is.
# Additionally, for modules that use GISTIC scores, we can just use the entire score
# Hence we will no longer run individual disease 

# The original OpenPBTA scripts were kept as is for references.

