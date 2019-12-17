
################## PURPOSE #################
# This script is to get rid of unwanted/trouble regions like telomeric, centromeric, and known seg-dups regions
# This script is to take in 2 files
#      1st: A reference/list of bad segments that contains the bad segments and their start/end positions
#      2nd: A list of CNVs that need the bad segments to be filtered out

################# ASSUMPTION ###############
# It is assumed that the reference file DO NOT have overlapping telomeric/centromeric/segments
# The provided reference file "bad_chromosomal_seg_updated_merged.txt" DOES NOT have
#   overlaping segments - thus the assumption is not violated.
# If the use of a different reference file is desired, make sure to merge the file first before using this script.

############################################

# Imports in the pep8 order https://www.python.org/dev/peps/pep-0008/#imports
# Standard library
import argparse
import subprocess
import sys
import os
import re

parser = argparse.ArgumentParser(description="""This script filters out telomeric,
                                             centromeric, and seg-dup regions""")
parser.add_argument('--reference', required=True,
                    help='path to the reference list of bad segments to be filtered out for')
parser.add_argument('--file', required=True,
                    help='path to the CNVs file that needs to be filtered')

args = parser.parse_args()

#querie_file = sys.argv[2]

#db_file = sys.argv[1]    #bad segment list


## Check to see if the CNVs file is empty, if so, make an empty variable
if os.stat(args.file).st_size == 0:
    final_content = []

## If the file has content, then load in the content
else:
    with open (args.file) as file:

        ## Strip the trailing spaces for each line
        stripped_content = [i.rstrip() for i in file.readlines()]

        ## Split each line up by any white space between columns
        final_content = [re.split('\s+',i) for i in stripped_content]

        ## If the 1st column has '1' instead of 'chr1' for the chromosome number, add in 'chr'
        if final_content[0][0].find('chr') == -1:
            for c,chromo in enumerate(final_content):
                final_content[c][0] = 'chr' + final_content[c][0]



## Works with the assumption that the segments in the file on the next line DO NOT OVERLAP
## Read in the reference file
with open(args.reference) as file:

    ## Strip the trailing spaces for each line
    reference_stripped = [i.rstrip() for i in file.readlines()]

    ## Split each line up by any white space between columns
    fin_reference_file = [re.split('\s+',i) for i in reference_stripped]

    ## If the 1st column has '1' instead of 'chr1' for the chromosome number, add in 'chr'
    if fin_reference_file[0][0].find('chr') == -1:
        for c,chromo in enumerate(fin_reference_file):
            fin_reference_file[c][0] = 'chr' + fin_reference_file[c][0]


## Define a variable to hold the final results
fin_list = []

## Loop through the CNVs file.
## If any segment in the reference file overlap >50%, then get rid of that CNV
for i in final_content:

    ## Initiate a variable for the total overlapping
    overlap = 0

    ## Get the start/end coordinates of the cnv
    start_cnv = int(i[1])
    end_cnv = int(i[2])

    ## Make sure end > start (end downstream from start)
    if start_cnv > end_cnv:
        start_cnv, end_cnv = end_cnv, start_cnv

    ## For each CNV, loop through the entire reference file
    for j in fin_reference_file:
        start_ref = int(j[1])
        end_ref = int(j[2])

        ## Make sure end > start
        if start_ref > end_ref:
            start_ref, end_ref = end_ref, start_ref

        ## If the chromosome number matches and the CNV segment overlaps with the reference segment
        if i[0] == j[0] and start_cnv <= end_ref and end_cnv >= start_ref:

            ## Get the overlapping coordinate
            start = max(start_cnv,start_ref)
            end = min(end_cnv,end_ref)

            ## Calculate the overlapping ratio
            ## If the CNV overlaps with many segments in the reference file
            ## We add the overlappting ratios -> This works due to the assumption on line 8
            overlap += (end - start + 1 )/(end_cnv-start_cnv + 1)


    ## If the overall overlapping ratio is smaller than 50%, we add the CNVs to a final list
    if overlap <= 0.5:
        fin_list = fin_list + [i]


## Print the list to STDOUT -> Python print() function prints to STDOUT by default.
try:
    for k in fin_list:
        print('\t'.join(k))
    sys.stderr.write('$$$ Filtering was successful\n')
except:
    sys.stderr.write('!!! Filtering was not successful\n')
