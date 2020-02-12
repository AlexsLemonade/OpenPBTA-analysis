#!/bin/bash

set -e 
set -o pipefail

# Configure environmental variables for MCR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/mcr/v83/runtime/glnxa64:/opt/mcr/v83/bin/glnxa64:/opt/mcr/v83/sys/os/glnxa64
export XAPPLRESDIR=/opt/mcr/v83/X11/app-defaults

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# The gzipped SEG file to use and a "nickname" for it 
# The FILEPREFIX is appended to the base file name for the
# gunzipped version of the file in scratch
SEGFILE=${SEGFILE:-"../../data/pbta-cnv-consensus.seg.gz"}
FILEPREFIX=${FILEPREFIX:-"consensus"}

# The results directory relative to where this script lives
# and the folder name that will contain the GISTIC results, respectively
# These will be constructed into a file path that is the basedir argument
# for GISTIC
RESULTSDIR=${RESULTSDIR:-"results"}
OUTPUTFOLDER=${OUTPUTFOLDER:-"pbta-cnv-consensus-gistic"}

basedir=${RESULTSDIR}/${OUTPUTFOLDER}
mkdir -p $basedir

# It's unclear if GISTIC can handle gzipped files
# So we will gunzip the file specified and put it in
# a "compartmentalized" scratch folder and name the file
# using the FILEPREFIX variable
scratchdir="../../scratch/uncompressed-seg-for-gistic"
mkdir -p $scratchdir
segfile=${scratchdir}/${FILEPREFIX}_seg_file_for_gistic.seg
gunzip -c $SEGFILE > $segfile

# This is the correct genome build for OpenPBTA
refgenefile=../../../gistic_install/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

# Run GISTIC! These parameters are from: https://github.com/d3b-center/OpenPBTA-workflows/blob/cb87a2b725d0e41d34a88436492830802c40f7f0/bash/run-gistic.sh#L16
../../../gistic_install/gp_gistic2_from_seg \
	-v 30 \
	-b $basedir \
	-seg $segfile \
	-refgene $refgenefile \
	-genegistic 1 \
	-smallmem 1 \
	-broad 1 \
	-twoside 1 \
	-brlen 0.98 \
	-conf 0.90 \
	-armpeel 1 \
	-savegene 1 \
	-gcm extreme \
	-js 2 \
	-rx 0
