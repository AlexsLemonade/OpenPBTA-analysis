#!/bin/bash

####################################
# This is script that reads the tcga-capture_kit-info.tsv and download
# all the uniq BED and add chr prefix to it.
#
# BAMs with `|` in the capture kit name and url are those with 
# more than one capture kit which neither the GDC nor its origin 
# data center could retrieve/figure out what the actual capture kit 
# had been applied. We should just generate intersec BED for those 
# samples and used that for our analysis.
#
# Those are all hg19 based coordinates, use UCSC online LiftOVer 
# to convert to GRCh38 based bed. Add chr to the TCGA bed to 
# make sure its format is liftover compatible.
####################################

TSV='results/tcga-capture_kit-info.tsv'

sed 1d $TSV \
| cut -f3 | tr "\|" "\n" | sort -u \
| while read i
do 
    filename=`basename $i`
    curl -s $i | awk '{print "chr"$0}' > ../../scratch/$filename
done
