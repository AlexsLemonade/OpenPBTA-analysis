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

Rscript --vanilla '01-snv-frequencies.R'

# Convert JSON to JSON Lines (JSONL) format with jq
#
# jq: https://stedolan.github.io/jq/ JSONL: https://jsonlines.org/
#
# Each 01-tpm-summary-stats.R output json file is an array of objects, i.e.
# [{k11:v11,k12:v12,...}, {k21:v21,k22:v22,...}, ...].
#
# jq docs:
#
# --compact-output / -c: By default, jq pretty-prints JSON output. Using this
# option will result in more compact output by instead **putting each JSON
# object on a single line**.
#
# Note putting each JSON object on a single line follows JSONL
#
# Array/Object Value Iterator: .[]
#
# If you use the .[index] syntax, but omit the index entirely, it will **return
# all of the elements of an array**. Running .[] with the input [1,2,3] will
# produce the numbers as three separate results, rather than as a single array.
#
# You can also use this on an object, and it will return all the values of the
# object.
#
# Commands adapted from
#
# - https://stackoverflow.com/a/48711608/4638182
# - https://stackoverflow.com/a/66709708/4638182
jq --compact-output '.[]' \
  results/snv-consensus-annotated-mut-freq.json \
  > results/snv-consensus-annotated-mut-freq.jsonl


# --no-name option stops the filename and timestamp from being stored in the
# output file. So rerun will have the same file.
gzip --no-name results/snv-consensus-annotated-mut-freq.json
gzip --no-name results/snv-consensus-annotated-mut-freq.jsonl
