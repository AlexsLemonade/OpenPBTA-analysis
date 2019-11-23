#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

## Run the python script to go from 1 big manta file, cnvkit file and freec file into 3 directories. Each directory with individual sample files. 

python3 src/scripts/merged_to_individual_files.py ../../data/pbta-sv-manta.tsv.gz ../../data/pbta-cnv-cnvkit.seg.gz ../../data/pbta-cnv-controlfreec.tsv.gz
