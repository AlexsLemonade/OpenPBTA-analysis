#!/bin/bash
# J. Taroni for CCDL 2019
# Create subset files for continuous integration

set -e
set -o pipefail

# Set defaults for release and biospecimen file name
BIOSPECIMEN_FILE=${BIOSPECIMEN_FILE:-biospecimen_ids_for_subset.RDS}
RELEASE=${RELEASE:-release-v5-20190924}

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

# get list of biospecimen ids for subset files
Rscript --vanilla 01-get_biospecimen_identifiers.R \
    --data_directory $FULL_DIRECTORY \
    --output_file $BIOSPECIMEN_FILE

# subset the files
Rscript --vanilla 02-subset_files.R \
  --biospecimen_file $BIOSPECIMEN_FILE \
  --output_directory $SUBSET_DIRECTORY

#### copy files that are not being subset --------------------------------------

# histologies file
cp $FULL_DIRECTORY/pbta-histologies.tsv $SUBSET_DIRECTORY

# all bed files
cp $FULL_DIRECTORY/*.bed $SUBSET_DIRECTORY

# if the md5sum.txt file already exists, get rid of it
rm -f $SUBSET_DIRECTORY/md5sum.txt
# create a new md5sum.txt file
cd $SUBSET_DIRECTORY
md5sum * > md5sum.txt

# Changelog does not get tracked
cd ../../..
cp $FULL_DIRECTORY/CHANGELOG.md $SUBSET_DIRECTORY
