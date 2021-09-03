#!/bin/bash
# Eric Wafula
# Pediatric OpenTargets 2021

# Purpose: Run a CNV consensus gene-level frequecies anaysis for PediatricOpenTargets

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# All file paths in this bash script are based on root directory of this Git repository,
# therefore the script should always be run from the root directory of OPenPedCan-analysis
# Adapted from OpenPBTA analysis modules

# Histology file path 
histology_file=data/histologies.tsv

# CNV consensus file path
cnv_file=data/consensus_wgs_plus_cnvkit_wxs.tsv.gz

# All cohorts independent primary tumor samples file path
all_cohorts_primary_tumors=data/independent-specimens.wgswxspanel.primary.tsv

# All cohorts independent relapse tumor samples file path
all_cohorts_relapse_tumors=data/independent-specimens.wgswxspanel.relapse.tsv

# Each cohort independent primary tumor samples file path
each_cohort_primary_tumors=data/independent-specimens.wgswxspanel.primary.eachcohort.tsv

# Each cohort independent relapse tumor samples file path
each_cohort_relapse_tumors=data/independent-specimens.wgswxspanel.relapse.eachcohort.tsv

####### compute CNV frequencies #############
python3 analyses/cnv-frequencies/01-cnv-frequencies.py $histology_file $cnv_file $all_cohorts_primary_tumors $all_cohorts_relapse_tumors $each_cohort_primary_tumors $each_cohort_relapse_tumors

####### compress the result files ######################
gzip analyses/cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq*

####### remove intermediate temporary and log files ######################
rm -f analyses/cnv-frequencies/results/annotator.log
rm -f analyses/cnv-frequencies/results/gene-level-cnv-consensus-mut-freq.tsv

