
"""
Create correlation matrices for polyA samples and ribodeplete (stranded) samples using gene expression data.
The correlations will be calculated using pairwise Spearman correlation of the RNA-Seq gene expression profiles.
The gene expression profile data will be filtered using a developed method described in Vaske et al.

https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/229
"""

import sys
import argparse
import pandas as pd
import utils


def run(input_matrix):
    # Load in data
    result = utils.read_rds(input_matrix)
    print("succesfully read in {} with {} rows".format(
        input_matrix, len(result)))

# TODO:
# run expression filter

# run variance filter

# do the rank sorting to make it Spearman

# run pearson correlation on sorted

# write output



def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-matrix", help="Path to gene expression matrix RDS file")
    args = p.parse_args()
    run(args.input_matrix)

if __name__ == "__main__":
    main()
