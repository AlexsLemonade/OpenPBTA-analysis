################## PURPOSE #################
# This script is to add column names/header to the consensus bed file

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
parser = argparse.ArgumentParser(description="""This script add columns names/headers to the consensus files""")
parser.add_argument('--file', required=True,
                    help='path to the bed file that needs the header')

args = parser.parse_args()


## Check to see if the bed file is empty, if so, make an empty dataframe
if os.stat(args.file).st_size == 0:
    new_mat = pd.DataFrame(columns=['a','b'])

## Use pandas to read in the bed file.
## Since Manta doesn't have copy numbers and "NA" is used for Manta's copy numbers, we need
## the "keep_default_na" so that the "NA"s don't get converted into NaN by pandas
else:
    new_mat = pd.read_csv(args.file, delimiter='\t', header=None, keep_default_na=False)

    ## Add a header to the DataFrame
    new_mat.columns = ['chrom','start','end','manta_CNVs','cnvkit_CNVs','freec_CNVs','CNV_type','sample','file_names']

## Output restructred file to stdout
new_mat.to_csv(sys.stdout, sep='\t', index=False, header=True)
