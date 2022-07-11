#!/bin/bash

# Jaclyn Taroni for ALSF CCDL 2020

# This script should always run as if it were being called from
# the directory it lives in.
cd "$(dirname "${BASH_SOURCE[0]}")"

# If this is CI, run the example included with GISTIC to check the installation
# The sample size for the subset files are too small otherwise
IS_CI=${OPENPBTA_CI:-0}

RESULTSDIR=results
DATADIR=../../data

# These will be constant for every disease
consensus_segfile="${RESULTSDIR}/pbta-cnv-consensus-gistic-only.seg.gz"
histologies_file="${DATADIR}/pbta-histologies.tsv"
filter_column="short_histology"

# Generate files that are compatible for GISTIC
# We should be able to check this step via CI
Rscript scripts/prepare_seg_for_gistic.R \
  --in_consensus ${DATADIR}/pbta-cnv-consensus.seg.gz \
  --out_consensus ${consensus_segfile} \
  --histology ${histologies_file}

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

  # run GISTIC for the whole cohort
  echo "Running GISTIC on the entire OpenPBTA cohort..."
  bash scripts/run-gistic-openpbta.sh

  # Now we'll run it on histologies with at least 100 WGS samples
  echo "Running GISTIC on specific histologies..."

  # directory where we will put the array list files
  # to a histology
  alf_directory="array_list_files"
  mkdir -p $alf_directory

  # These will be constant for every disease
  filter_column="short_histology"

  # borrowed from Josh Shapiro in analyses/interaction plot
  # associative array of diseases to test; chosen by those that are most common
  # in the openPBTA dataset
  # all of these histologies have >100 WGS samples
  declare -A disease
  disease[LGAT]="LGAT"
  disease[HGAT]="HGAT"
  disease[Medulloblastoma]="Medulloblastoma"

  for disease_id in "${!disease[@]}"; do
    echo "    $disease_id"

    # generate a subset SEG file for this disease
    array_list_file="${alf_directory}/${disease_id,,}-array-list.txt"
    Rscript --vanilla scripts/generate-array-file.R \
      --segfile $consensus_segfile \
      --metadata $histologies_file \
      --filter_column $filter_column \
      --filter_value ${disease[$disease_id]} \
      --output_file $array_list_file

    ARRAYLIST=../${array_list_file} \
    FILEPREFIX=$disease_id \
    OUTPUTFOLDER=pbta-cnv-consensus-${disease_id,,}-gistic \
    bash scripts/run-gistic-openpbta.sh
   done
fi
