#!/bin/bash
set -eo pipefail

# set the analyses folder as workdir variable
#WORKDIR='analyses/tcga-capture-kit-investigation'
WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# retrieve the TCGA's exome capture kit from GDC's file API endpoint
# this will generate tcga-capture_kit-info.tsv file under the result folder
python3 $WORKDIR/scripts/get-tcga-capture_kit.py

# download all BED from tcga-capture_kit-info.tsv and add chr prefix
sed 1d $WORKDIR/results/tcga-capture_kit-info.tsv \
| cut -f3 | tr "\|" "\n" | sort -u \
| while read i
do 
    filename=`basename $i`
    curl -s $i | awk '{print "chr"$0}' > ../../scratch/$filename
done

##Downloading chain.gz file for CrossMap command 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz  -O ../../scratch/hg19ToHg38.over.chain.gz

## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
for i in `ls ../../scratch/*.targetIntervals.bed`; do out=$(echo $i | sed 's/.bed/.Crossmapped_to_Gh38.bed/g'); CrossMap.py bed ../../scratch/hg19ToHg38.over.chain.gz $i $out; done

## This command sorts and merges all BED files in case there  are any overlaps in the BED files
for i in `ls ../../scratch/*.targetIntervals.Crossmapped_to_Gh38.bed`; do out=$(echo $i | sed 's/Crossmapped_to_Gh38/Gh38/g'); bedtools sort -i $i | bedtools merge -i - > $out;mv $out $WORKDIR/results; done
rm ../../scratch/*.targetIntervals.Crossmapped_to_Gh38.bed

## get intersection between all the BED files
bedtools intersect -a ../../data/WGS.hg38.strelka2.unpadded.bed -b ../../data/WGS.hg38.mutect2.vardict.unpadded.bed \
| bedtools intersect -a - -b ../../data/WGS.hg38.lancet.unpadded.bed \
| bedtools intersect -a - -b $WORKDIR/results/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.Gh38.bed \
| bedtools intersect -a - -b $WORKDIR/results/whole_exome_agilent_designed_120.targetIntervals.Gh38.bed \
| bedtools intersect -a - -b $WORKDIR/results/whole_exome_agilent_plus_tcga_6k.targetIntervals.Gh38.bed \
| bedtools intersect -a - -b $WORKDIR/results/tcga_6k_genes.targetIntervals.Gh38.bed \
> ../../scratch/intersect-all-tcga_pbta.bed
