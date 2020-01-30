#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


SCRATCHDIR=../../scratch/copy_consensus
# make directories:
mkdir -p $SCRATCHDIR
mkdir -p results

## Run the python script to go from 1 big manta file, cnvkit file and freec file into 3 directories. 
## Each directory with individual sample files.

python3 scripts/merged_to_individual_files.py \
    --manta ../../data/pbta-sv-manta.tsv.gz \
    --cnvkit ../../data/pbta-cnv-cnvkit.seg.gz \
    --freec ../../data/pbta-cnv-controlfreec.tsv.gz \
    --snake $SCRATCHDIR/config_snakemake.yaml \
    --scratch $SCRATCHDIR \
    --uncalled results/uncalled_samples.tsv


## Run the Snakemake pipeline
## This Snakemake is to produce copy number consensus, it:
## 1) Filters out the CNVs results of the 3 call methods
## 2) Performs a reciprocal comparison between 2 call methods to find common CNVs agreed upon by those 2 methods
## 3) Repeats step 2 for all pairs made from the 3 call methods 
## 4) Finally merges the consensus calls together into one big consensus file

## The snakemake flag options are:
## -s : Point to the location of the Snakemake file
## --configfile : Point to the location of the config file
## -j : Set available cores, in this case, when no number is provided, thus use all available cores
## -p : Print shell command that will be executed
## --restart-times : Define the times a job restarts when run into an error before giving up
## --latency-wait: Define the number of seconds to wait for a file to show up after that file has been created

snakemake \
    -s Snakefile \
    --configfile $SCRATCHDIR/config_snakemake.yaml \
    -j \
    --restart-times 2

