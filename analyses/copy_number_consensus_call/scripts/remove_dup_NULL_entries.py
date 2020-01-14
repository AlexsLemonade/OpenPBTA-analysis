################## PURPOSE #################
# This script is to remove duplicated coordinates reulted from Bedtools merging
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


def remove_dup_and_null(cell):
    """ Remove duplicated coordinate sets and NULLs from a comma separated string """

    ## Split the cell by commas and convert it to a set, this removes duplicated elements
    set_cell = set(cell.split(','))

    ## Keep NULLs if it's the only entry in a cell, else, remove NULLs
    if "NULL" in set_cell and len(set_cell) > 1:
        set_cell.remove("NULL")

    ## Join the set
    if len(set_cell) > 1:
        sorted_cell = sorted(set_cell, key=lambda x: int(x.split(":")[0]))
    else:
        sorted_cell = set_cell
    new_cell = ",".join(sorted_cell)
    return new_cell


## Define parser for the input file
parser = argparse.ArgumentParser(
    description="""This script removes duplicates coordinate sets and NULLs from the individual callers columns"""
)
parser.add_argument('--file', required=True,
                    help='path to the file that needs duplicates and NULLs removed')

args = parser.parse_args()

## Read in the input file
new_mat = pd.read_csv(args.file, delimiter='\t', header=[0], keep_default_na=False)

## Take out the three columns to be operated on and apply the "remove_dup_and_null" function on those columns
new_mat[['manta_CNVs','cnvkit_CNVs','freec_CNVs']] = \
    new_mat[['manta_CNVs','cnvkit_CNVs','freec_CNVs']].applymap(remove_dup_and_null)


## Output restructred file to stdout
new_mat.to_csv(sys.stdout, sep='\t', index=False, header=True)
