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
autosomes_cnv_file=data/consensus_wgs_plus_cnvkit_wxs_autosomes.tsv.gz
allosomes_cnv_file=data/consensus_wgs_plus_cnvkit_wxs_x_and_y.tsv.gz

# Independent primary tumor samples file path
primary_tumors=analyses/independent-samples/results/independent-specimens.wgswxspanel.primary.eachcohort.tsv

# Independent relapse tumor samples file path
relapse_tumors=analyses/independent-samples/results/independent-specimens.wgswxspanel.relapse.eachcohort.tsv

####### compute autosomes CNV frequencies #############
python3 analyses/cnv-frequencies/01-cnv-frequencies.py $histology_file $autosomes_cnv_file $primary_tumors $relapse_tumors

####### compute allosomes CNV frequencies #############
python3 analyses/cnv-frequencies/01-cnv-frequencies.py $histology_file $allosomes_cnv_file $primary_tumors $relapse_tumors

####### merge result files ######################
# intricate bash command concatenates the autosomes and allosomes annotated CNV frequencies TSV files, and excludes the header line for the second 
# file to avoid duplicate header in the merged file
# bash command "tail -n +2" outputs content of file starting from the second line
# bash Process substitution operator "<(...)" treats the sequence of enclosed commands as file 
cat analyses/cnv-frequencies/results/consensus_wgs_plus_cnvkit_wxs_autosomes_annot_freq.tsv <(tail -n +2 analyses/cnv-frequencies/results/consensus_wgs_plus_cnvkit_wxs_x_and_y_annot_freq.tsv)> analyses/cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq.tsv
# simple bash command concatenates the autosomes and allosomes annotated CNV frequencies JSONL files
cat analyses/cnv-frequencies/results/consensus_wgs_plus_cnvkit_wxs_autosomes_annot_freq.jsonl analyses/cnv-frequencies/results/consensus_wgs_plus_cnvkit_wxs_x_and_y_annot_freq.jsonl > analyses/cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq.jsonl

####### compress the final merged result files ######################
gzip analyses/cnv-frequencies/results/gene-level-cnv-consensus-annotated-mut-freq*

####### remove intermediate temporary and log files ######################
rm -f analyses/cnv-frequencies/results/consensus_wgs_plus_cnvkit_wxs_*
rm -f analyses/cnv-frequencies/results/annotator.log

