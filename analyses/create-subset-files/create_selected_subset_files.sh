#!/bin/bash
# J. Taroni for CCDL 2019
# Create subset files for CI for selected files in a new release using the
# old release files and biospecimen ID lists (stored as RDS).
# Use this script if you have the old release subset files available in
# ../../data/testing/<old release>

set -e
set -o pipefail

# Set defaults for release and biospecimen file name
BIOSPECIMEN_FILE=${BIOSPECIMEN_FILE:-biospecimen_ids_for_subset.RDS}
OLD_RELEASE=${OLD_RELEASE:-release-v7-20191031}
NEW_RELEASE=${NEW_RELEASE:-release-v8-20191104}
SELECTED_FILES=${SELECTED_FILES:-NULL}

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

OLD_DIRECTORY=../../data/testing/${OLD_RELEASE}
NEW_DIRECTORY=../../data/testing/${NEW_RELEASE}
FULL_DIRECTORY=../../data/${NEW_RELEASE}

# copy all the old files from the old testing directory to the new testing
# directory
cp $OLD_DIRECTORY/* $NEW_DIRECTORY

# now run the subsetting on the selected files
# subset the files
Rscript --vanilla 02-subset_files.R \
  --biospecimen_file $BIOSPECIMEN_FILE \
  --output_directory $NEW_DIRECTORY \
  --selected_files $SELECTED_FILES \
  --new_release $NEW_RELEASE

# histologies file
\cp $FULL_DIRECTORY/pbta-histologies.tsv $NEW_DIRECTORY

# independent specimen files
\cp $FULL_DIRECTORY/independent-specimens*.tsv $NEW_DIRECTORY

# all bed files
\cp $FULL_DIRECTORY/*.bed $NEW_DIRECTORY

# if the md5sum.txt file already exists, get rid of it
rm -f $NEW_DIRECTORY/md5sum.txt
# create a new md5sum.txt file
cd $NEW_DIRECTORY
md5sum * > md5sum.txt

# Changelog does not get tracked
cd "$script_directory"
\cp $FULL_DIRECTORY/CHANGELOG.md $NEW_DIRECTORY
