#!/bin/bash

# Jaclyn Taroni for ALSF CCDL 2020

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# If this is CI, run the example included with GISTIC
# The sample size for the subset files are too small otherwise
IS_CI=${OPENPBTA_CI:-0}

if [[ "$IS_CI" -gt "0" ]]
then
  cd /home/rstudio/gistic_install && ./run_gistic_example
else

	# run GISTIC for the whole cohort 
	echo "Running GISTIC on the entire OpenPBTA cohort..."
	bash scripts/run-gistic-openpbta.sh

	# Now we'll run it on histologies with at least 100 WGS samples
	echo "Running GISTIC on specific histologies..."

	# directory where we will put the SEG files that are specific 
	# to a histology
	subset_directory="seg_files"
	mkdir -p $subset_directory

	# These will be constant for every disease
	consensus_segfile="../../data/pbta-cnv-consensus.seg.gz"
	histologies_file="../../data/pbta-histologies.tsv"
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
		subset_seg_file="${subset_directory}/${disease_id,,}_consensus_seg.gz"
		Rscript --vanilla scripts/subset-seg-file.R \
			--segfile $consensus_segfile \
			--metadata $histologies_file \
			--filter_column $filter_column \
			--filter_value $disease_id \
			--output_file $subset_seg_file

		SEGFILE=../${subset_seg_file} \
		FILEPREFIX=$disease_id \
		OUTPUTFOLDER=pbta-cnv-consensus-${disease_id,,}-gistic \
		bash scripts/run-gistic-openpbta.sh

	done

fi 