#!/bin/bash


## Run the python script to go from 1 big manta file, cnvkit file and freec file into 3 directories. Each directory with individual sample files. 

python3 src/scripts/merged_to_individual_files.py ../../data/pbta-sv-manta.tsv.gz ../../data/pbta-cnv-cnvkit.seg.gz ../../data/pbta-cnv-controlfreec.tsv.gz
