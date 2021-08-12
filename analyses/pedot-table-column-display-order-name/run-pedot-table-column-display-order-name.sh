#!/bin/bash
# PediatricOpenTargets 2021
# Yuanchao Zhang
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
# copied from the run_in_ci.sh file at
# <https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/scripts/>
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

mkdir -p results

echo 'Run 01-generate-pedot-column-display-order-name-xlsx.py...'
python3 01-generate-pedot-column-display-order-name-xlsx.py

echo 'Done running 01-generate-pedot-column-display-order-name-xlsx.py.'

echo 'Done running run-pedot-table-column-display-order-name.sh'
