################## PURPOSE #################
# This script removes duplicated coordinates that resulted from Bedtools merging
# and removes NULLs if they are not stand-alone values in a cell

################# ASSUMPTION ###############


############################################

# Imports in the pep8 order https://www.python.org/dev/peps/pep-0008/#imports
# Standard library
import argparse
import sys

# Related third party
import numpy as np
import pandas as pd


def cnv_overlap(cnv, start, end):
    #check if a cnv overlaps the interval defined by start and end
    if cnv == "NULL":
        return(True)
    cnv_start, cnv_end, _ = cnv.split(':', maxsplit = 2)
    return(int(cnv_end) > start and int(cnv_start) < end)

def remove_dup_null_outside(cell, start, end):
    """ Remove duplicated coordinate sets and NULLs from a comma separated string """

    ## Split the cell by commas and convert it to a set, this removes duplicated elements
    cnv_set = set(cell.split(','))

    ## Keep NULLs if it's the only entry in a cell, else, remove NULLs
    if "NULL" in cnv_set and len(cnv_set) > 1:
        cnv_set.remove("NULL")

    ## remove cells that are not in the interval (can happen when resolving overlaps)
    cnv_set = list(filter(lambda cnv: cnv_overlap(cnv, start, end), cnv_set))

    ## Join the set
    if len(cnv_set) > 1:
        sorted_cnvs = sorted(cnv_set, key=lambda x: int(x.split(":")[0]))
    else:
        sorted_cnvs = cnv_set
    new_cell = ",".join(sorted_cnvs)
    return new_cell


## Define parser for the input file
parser = argparse.ArgumentParser(
    description="""This script removes duplicate coordinate sets and NULLs from the individual callers columns"""
)
parser.add_argument('--file', required=True,
                    help='path to the file that needs duplicates and NULLs removed')

args = parser.parse_args()

## Read in the input file
new_mat = pd.read_csv(args.file, delimiter='\t', header=[0], keep_default_na=False)

## Take out the three columns to be operated on and apply the "remove_dup_and_null" function on those columns
new_mat['manta_CNVs'] = \
    [remove_dup_null_outside(c, s, e) for c, s, e in zip(new_mat['manta_CNVs'], new_mat['start'], new_mat['end'])]

new_mat['cnvkit_CNVs'] = \
    [remove_dup_null_outside(c, s, e) for c, s, e in zip(new_mat['cnvkit_CNVs'], new_mat['start'], new_mat['end'])]

new_mat['freec_CNVs'] = \
    [remove_dup_null_outside(c, s, e) for c, s, e in zip(new_mat['freec_CNVs'], new_mat['start'], new_mat['end'])]


## Output restructred file to stdout
new_mat.to_csv(sys.stdout, sep='\t', index=False, header=True)
