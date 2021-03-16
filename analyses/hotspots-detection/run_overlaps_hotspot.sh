#!/bin/bash

# K S Gaonkar
set -e
set -o pipefail



# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

cancer_hotspot_file=input/hotspots_v2.xls
genomic_site_hotspot_file=input/tert_promoter_hotspots.tsv

# find hotspot overlaps in strelka2
Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-strelka2.vep.maf.gz \
        --caller strelka2 \
        --cancer_hotspot_file $cancer_hotspot_file \
        --genomic_site_hotspot_file $genomic_site_hotspot_file 

# find hotspot overlaps in mutect2
Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-mutect2.vep.maf.gz \
        --caller mutect2 \
        --cancer_hotspot_file $cancer_hotspot_file \
        --genomic_site_hotspot_file $genomic_site_hotspot_file

# find hotspot overlaps in lancet
Rscript 00-subset-maf.R --maffile ../../data/pbta-snv-lancet.vep.maf.gz \
        --caller lancet \
        --cancer_hotspot_file $cancer_hotspot_file \
        --genomic_site_hotspot_file $genomic_site_hotspot_file


# find hotspot overlaps in vardict
tmpdir=../../scratch/hotspots-detection
mkdir -p $tmpdir
gunzip -c ../../data/release-v18-20201123/pbta-snv-vardict.vep.maf.gz | split -l 10000000 - ${tmpdir}/pbta-snv-vardict.vep.maf.

n=1
for file in $(ls ${tmpdir}/pbta-snv-vardict.vep.maf.*); do
 gzip $file ; 
 Rscript 00-subset-maf.R --maffile ${file}.gz\
        --caller vardict0${n} \
        --cancer_hotspot_file $cancer_hotspot_file \
        --genomic_site_hotspot_file $genomic_site_hotspot_file ;
 n=$(( $n + 1 ));
done
