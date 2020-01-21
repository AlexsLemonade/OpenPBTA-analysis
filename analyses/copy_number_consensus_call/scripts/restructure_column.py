################## PURPOSE #################
# This script is to restructure the 4th(start pos), 5th(end pos), 6th(copy number), and 7th(seg.mean) columns of the input file
# Because these columns hold the "raw" coordinates of the CNV after filtering out for IGLL, telomerics, centromeric
# and segup regions, we want to retain these info for any analysis further down stream.
# Thus this script format columns 4, 5, 6, and 7 of a CNV into ---  start_pos:end_pos:copy_number:seg_mean  ----
# So if a CNV is part of a FINAL consensus CNV, its (start_pos:end_pos:copy_number:seg_mean) information will
# be carried through and contained in the FINAL FILE

# For example
#          column 4                   column 5             column 6        column 7
# 2777350,6035598,6807639       6035597,6807638,7659545     3,3,3         0.5,0.2,0.4

# Will become 2777350:6035597:3:0.5,6035598:6807638:3:0.2,6807639:7659545:3:0.4,
################# ASSUMPTION ###############

# The bed file has no header

############################################

# Imports in the pep8 order https://www.python.org/dev/peps/pep-0008/#imports
# Standard library
import argparse
import sys
import os

# Related third party
import numpy as np
import pandas as pd

## Define parser for the input file
parser = argparse.ArgumentParser(description="""This script restructure the 4th, 5th, and 6th columns""")
parser.add_argument('--file', required=True,
                    help='path to the bed file that needs to be restructred')

args = parser.parse_args()


## Check to see if the bed file is empty, if so, make an empty dataframe
if os.stat(args.file).st_size == 0:
    new_mat = pd.DataFrame(columns=['a','b'])

## Use pandas to read in the bed file.
## Since Manta doesn't have copy numbers and "NA" is used for Manta's copy numbers, we need
## the "keep_default_na" so that the "NA"s don't get converted into NaN by pandas
else:
    file = pd.read_csv(args.file, delimiter='\t', header=None, keep_default_na=False)

    ## Make a new column to store info in the new format (start_pos:end_pos:copy_number)
    new_column = file.iloc[:,(0)].copy()

    ## For every CNV in the file, reformat columns 4, 5, and 6
    for i_index,i in enumerate(file.itertuples()):
        ## Since bed file uses a comma as default for collapsing information, we split the information using ','
        list_start = str(i._4).split(',')
        list_end = str(i._5).split(',')
        list_cn = str(i._6).split(',')
        list_segmean = str(i._7).split(',')

        ## Rearange the split information into the new format
        cnv_list = ''.join(['{}:{}:{}:{},'.format(*cnv) for cnv in zip(list_start, list_end, list_cn, list_segmean)])


        ## Add the list to the new column
        new_column.iloc[i_index] = cnv_list


    ## Make a new file by combining columns 1,2,3, the new column, and the copy number type column.
    new_mat = (pd.concat([file.iloc[:,0:3].copy(),new_column,file.iloc[:,-1].copy()], axis=1))

## Output restructred file to stdout
new_mat.to_csv(sys.stdout, sep='\t', index=False, header=False)
