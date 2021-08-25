#!/bin/bash
# J. Taroni for CCDL 2019
# Create subset files for continuous integration

set -e
set -o pipefail

# Set defaults for release and biospecimen file name
BIOSPECIMEN_FILE=${BIOSPECIMEN_FILE:-biospecimen_ids_for_subset.RDS}
RELEASE=${RELEASE:-release-v21-20210820}
NUM_MATCHED=${NUM_MATCHED:-15}

# This option controls whether or not the two larger MAF files are skipped as
# part of subsetting -- the idea is that setting RUN_LOCAL=1 will allow for
# local testing
RUN_LOCAL=${RUN_LOCAL:-0}

# Use SKIP_SUBSETTING=1 to skip the subsetting steps and only copy full files
# and generate a new md5sum.txt file - this can be useful if the only files
# getting updated in a release are those that are copied in full
SKIP_SUBSETTING=${SKIP_SUBSETTING:-0}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# directories that hold the full files for the release and the subset files
# generated via these scripts
FULL_DIRECTORY=../../data/$RELEASE
SUBSET_DIRECTORY=../../data/testing/$RELEASE

#### generate subset files -----------------------------------------------------

if [ "$SKIP_SUBSETTING" -lt "1" ]; then

  # get list of biospecimen ids for subset files
  Rscript --vanilla 01-get_biospecimen_identifiers.R \
      --data_directory $FULL_DIRECTORY \
      --output_file $BIOSPECIMEN_FILE \
      --num_matched $NUM_MATCHED \
      --local $RUN_LOCAL

  # subset the files
  Rscript --vanilla 02-subset_files.R \
    --biospecimen_file $BIOSPECIMEN_FILE \
    --output_directory $SUBSET_DIRECTORY

fi

#### copy files that are not being subset --------------------------------------

# histologies file
cp $FULL_DIRECTORY/pbta-histologies.tsv $SUBSET_DIRECTORY

# base histologies file
cp $FULL_DIRECTORY/pbta-histologies-base.tsv $SUBSET_DIRECTORY

# recurrently fused genes by histologies file
cp $FULL_DIRECTORY/pbta-fusion-recurrently-fused-genes-byhistology.tsv $SUBSET_DIRECTORY

# GISTIC output
cp $FULL_DIRECTORY/pbta-cnv-cnvkit-gistic.zip $SUBSET_DIRECTORY
cp $FULL_DIRECTORY/pbta-cnv-consensus-gistic.zip $SUBSET_DIRECTORY

# independent specimen files
cp $FULL_DIRECTORY/independent-specimens*.tsv $SUBSET_DIRECTORY

# all bed files
cp $FULL_DIRECTORY/*.bed $SUBSET_DIRECTORY

# data file description
cp $FULL_DIRECTORY/data-files-description.md $SUBSET_DIRECTORY

# TCGA related files
cp $FULL_DIRECTORY/pbta-tcga* $SUBSET_DIRECTORY
cp $FULL_DIRECTORY/tcga-snv* $SUBSET_DIRECTORY

# STAR logs
cp $FULL_DIRECTORY/pbta-star* $SUBSET_DIRECTORY

# MEND QC files
cp $FULL_DIRECTORY/pbta-mend* $SUBSET_DIRECTORY

# fusion summary files
cp $FULL_DIRECTORY/fusion_summary* $SUBSET_DIRECTORY

# MB pathology subtypes
cp $FULL_DIRECTORY/pbta-mb-pathology-subtypes.tsv $SUBSET_DIRECTORY

# if the md5sum.txt file already exists, get rid of it
cd $SUBSET_DIRECTORY
rm -f md5sum.txt
# create a new md5sum.txt file
md5sum * > md5sum.txt
cd "$script_directory" || exit

# the release notes and not included in md5sum.txt
cp $FULL_DIRECTORY/release-notes.md $SUBSET_DIRECTORY
