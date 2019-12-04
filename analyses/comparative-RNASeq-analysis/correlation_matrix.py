
"""
Create correlation matrices for polyA samples and ribodeplete (stranded) samples using gene expression data.
The correlations will be calculated using pairwise Spearman correlation of the RNA-Seq gene expression profiles.
The gene expression profile data will be filtered using a developed method described in Vaske et al.

https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/229
"""

import sys
import argparse
import math
import numpy as np
import pandas as pd
import scipy.stats
import sklearn.metrics.pairwise as sklp
import utils

def prepare_expression(input_matrix,
    proportion_unexpressed=0.8,
    variance_filter_level=0.2,
    verbose=False,
    prefix=""):
    """Loads samples, converts to log2(TPM+1), and applies expression and variance filters.
    returns dataframe of converted and filtered expression.
    writes file: scratch/{prefix}filtered_genes_to_keep.rds"""

    print_v = print if verbose else lambda *a, **k: None
    filtered_genelist_filename = "{}filtered_genes_to_keep.rds".format(prefix)

    print_v("Loading samples {}".format(input_matrix))
    raw_samples = utils.read_rds(input_matrix, location="data")

    # Convert to log2(tpm+1)
    samples = raw_samples.apply(lambda x: np.log2(x+1))

    # run expression filter (filter_out_genes_unexpressed_in_most_samples.py)
    # Remove any genes that have 0 expression in more samples than proportion_unexpressed
    print_v("Running expression filter")
    max_ok_zeroes = len(samples.columns) * proportion_unexpressed
    has_few_enough_zeroes = samples.apply(
        lambda s: len(s[s<=0]) < max_ok_zeroes,
        axis=1)
    expression_filtered_df = samples[has_few_enough_zeroes]

    # run variance filter (filter_out_lowest_varying_genes.py)
    # Sort remaining genes by variance and drop the least-varying
    # How many to drop is controlled by variance-filter_level
    print_v("Running variance filter")
    variance = expression_filtered_df.apply(np.std, axis=1)
    cut_proportion = int(math.ceil(len(variance)*variance_filter_level))
    keep_proportion = len(variance) - cut_proportion

    expr_var_filtered_genelist = pd.DataFrame(variance.nlargest(keep_proportion).index)
    expr_var_filtered_df = expression_filtered_df.filter(expr_var_filtered_genelist["gene_id"], axis=0)

    utils.write_rds(expr_var_filtered_genelist, filtered_genelist_filename, location="scratch")

    return expr_var_filtered_df


def calculate_correlation(expr_var_filtered_df, verbose=False, prefix=""):
    print_v = print if verbose else lambda *a, **k: None
    """Takes an expression matrix dataframe and calculates spearman correlations.
    output file: results/{prefix}all_by_all_correlations.tsv"""

    all_by_all_filename = "{}all_by_all_correlations.rds".format(prefix)

    # spearman is a rank transformed pearson ('correlation') metric
    # Rank transform to get ordering of genes within sample
    print_v("Rank transforming")
    filtered_rank_transformed_df = np.apply_along_axis(scipy.stats.rankdata,0,expr_var_filtered_df)

    # run pearson correlation on ranked genes. Note that df must have samples=rows.
    print_v("Running pairwise distances")
    x_corr = 1 - sklp.pairwise_distances(X=filtered_rank_transformed_df.transpose(), metric="correlation")

    # Reapply the axis labels
    labels=expr_var_filtered_df.columns
    all_by_all_df = pd.DataFrame(x_corr, index=labels,columns=labels)

    print_v("Writing to file")
    utils.write_rds(all_by_all_df, all_by_all_filename, location="results")
    return all_by_all_df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-matrix", help="Path to gene expression matrix RDS file")
    p.add_argument("--prefix", help="Prefix for output filenames")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()
    expression = prepare_expression(args.input_matrix, verbose=args.verbose, prefix=args.prefix)
    calculate_correlation(expression, verbose=args.verbose, prefix=args.prefix)


if __name__ == "__main__":
    main()
