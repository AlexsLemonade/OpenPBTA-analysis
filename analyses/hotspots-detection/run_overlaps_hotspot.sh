#!/bin/bash

# K S Gaonkar
set -e
set -o pipefail



# This script should always run as if it were being called from
# the directory it lives in.
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

cancer_hotspot_folder=input/hotspots_database
genomic_site_hotspot_file=input/tert_promoter_hotspots.tsv

# find hotspot overlaps in strelka2
Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-strelka2.vep.maf.gz \
        --caller strelka2 \
        --cancer_hotspot_folder $cancer_hotspot_folder \
        --genomic_site_hotspot_file $genomic_site_hotspot_file 

