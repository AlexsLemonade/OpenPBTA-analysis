#!/bin/bash

# K S Gaonkar
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# input hotspots folder and files
cancer_hotspot_folder=input/hotspots_database
genomic_site_hotspot_file=input/tert_promoter_hotspots.tsv

# Create the maf_coltypes.RDS 
# this will help in reading maf files correctly in following scripts
Rscript maf_coltypes.R

# tmp directory in scratch folder to save intermediate filtered maf files
tmp_dir="../../scratch/hotspot-detection"

mkdir -p $tmp_dir

# given a list of callers 
for caller in strelka2 mutect2 lancet vardict; do
    # create a new file to store gene filtered maf 
    rm -f $tmp_dir/pbta-snv-${caller}.vep.genefiltered.maf
    touch $tmp_dir/pbta-snv-${caller}.vep.genefiltered.maf
    # create unique list of genes in hotspot files to filter
    echo "Hugo_Symbol" > $tmp_dir/hotspot_genes.txt
    cat ${cancer_hotspot_folder}/hotspot_database_2017_snv.tsv ${cancer_hotspot_folder}/hotspot_database_2017_indel.tsv|cut -f 1| grep -v "Hugo_Symbol"|sort|uniq >> $tmp_dir/hotspot_genes.txt
    # filter maf for genes in mskcc hotspot tsv files
    # if FNR==NR which will only be the case while reading $tmp_dir/hotspot_genes.txt
    # make an array `arr` with the gene names
    # if $1 in ../../data/pbta-snv-${caller}.vep.maf.gz in array of genes `arr` print the row
    awk 'FNR==NR {arr[$1];next} $1 in arr' $tmp_dir/hotspot_genes.txt <(gunzip -c ../../data/pbta-snv-${caller}.vep.maf.gz)  >> $tmp_dir/pbta-snv-${caller}.vep.genefiltered.maf
    # filter for TERT for non-coding hotspot
    gunzip -c ../../data/pbta-snv-${caller}.vep.maf.gz|awk -F "\t" '{ if($1=="TERT") print}' >> $tmp_dir/pbta-snv-${caller}.vep.genefiltered.maf
    # find hotspot overlaps
    Rscript 00-subset-maf.R --maffile ${tmp_dir}/pbta-snv-${caller}.vep.genefiltered.maf \
        --caller ${caller} \
        --cancer_hotspot_folder $cancer_hotspot_folder \
        --genomic_site_hotspot_file $genomic_site_hotspot_file \
        --output_file  ${tmp_dir}/${caller}_hotspots.RDS
done

# combine hotspots with consensus 
Rscript -e "rmarkdown::render('01-create-hotspot-maf.Rmd')"

