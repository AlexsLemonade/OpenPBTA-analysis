
"""
Create correlation matrices for polyA samples and ribodeplete (stranded) samples using gene expression data.
The correlations will be calculated using pairwise Spearman correlation of the RNA-Seq gene expression profiles.
The gene expression profile data will be filtered using a developed method described in Vaske et al.

https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/229
"""

import os
import sys
import argparse
import math
import numpy as np
import pandas as pd
import scipy.stats
import sklearn.metrics.pairwise as sklp
import utils

def prepare_expression(input_matrix,
                       scratch_dir,
                       normalized_samples_path,
                       proportion_unexpressed=0.8,
                       variance_filter_level=0.2,
                       verbose=False,
                       prefix=""):
    """Loads samples, converts to log2(TPM+1), and applies expression and variance filters.
    returns dataframe of converted and filtered expression.
    input_matrix = path to .rds file
    writes file: scratch/{prefix}filtered_genes_to_keep.rds"""

    print_v = print if verbose else lambda *a, **k: None
    filtered_genelist_filepath = os.path.join(scratch_dir, "{}filtered_genes_to_keep.rds".format(prefix))

    print_v("Loading samples {}".format(input_matrix))
    raw_samples = utils.read_rds(input_matrix)

    # Convert to log2(tpm+1)
    samples = np.log2(raw_samples+1)
    print_v("Writing normalized samples to {}".format(normalized_samples_path))
    utils.write_rds(samples, normalized_samples_path)

    # run expression filter
    # Remove any genes that have 0 expression in more samples than proportion_unexpressed
    print_v("Running expression filter")
    max_ok_zeroes = len(samples.columns) * proportion_unexpressed
    has_few_enough_zeroes = samples.apply(
        lambda s: len(s[s<=0]) < max_ok_zeroes,
        axis=1)
    expression_filtered_df = samples[has_few_enough_zeroes]

    # run variance filter
    # Sort remaining genes by variance and drop the least-varying
    # How many to drop is controlled by variance_filter_level
    print_v("Running variance filter")
    variance = expression_filtered_df.apply(np.std, axis=1)
    cut_proportion = int(math.ceil(len(variance)*variance_filter_level))
    keep_proportion = len(variance) - cut_proportion

    expr_var_filtered_genelist = pd.DataFrame(variance.nlargest(keep_proportion).index)
    expr_var_filtered_df = expression_filtered_df.filter(expr_var_filtered_genelist["gene_id"], axis=0)

    utils.write_rds(expr_var_filtered_genelist, filtered_genelist_filepath)

    return expr_var_filtered_df


def calculate_correlation(expr_var_filtered_df, scratch_dir, verbose=False, prefix=""):
    """Takes an expression matrix dataframe and calculates spearman correlations.
    (spearman is a rank transformed pearson ('correlation') metric)
    output file: scratch/{prefix}all_by_all_correlations.tsv"""

    print_v = print if verbose else lambda *a, **k: None

    all_by_all_filepath = os.path.join(scratch_dir, "{}all_by_all_correlations.rds".format(prefix))

    # Rank transform - values are now that rank of that gene within each sample
    print_v("Rank transforming")
    filtered_rank_transformed_df = np.apply_along_axis(scipy.stats.rankdata, 0, expr_var_filtered_df)

    # Run pearson metric on ranked genes. For pairwise_distances we need samples=rows so transpose.
    print_v("Running pairwise distances")
    x_corr = 1 - sklp.pairwise_distances(X=filtered_rank_transformed_df.transpose(), metric="correlation")

    # Reapply the axis labels
    labels=expr_var_filtered_df.columns
    all_by_all_df = pd.DataFrame(x_corr, index=labels, columns=labels)

    print_v("Writing to file {}".format(all_by_all_filepath))
    utils.write_rds(all_by_all_df, all_by_all_filepath)
    return all_by_all_df

def main():
    """Creates correlation matrices for dataset."""
    # This script should always run as if it were being called from
    # the directory it lives in.
    os.chdir(sys.path[0])

    p = argparse.ArgumentParser()
    p.add_argument("input_path", metavar="input-path", help="Path to input Rds file.")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--scratch", 
                   default=os.path.join("..", "..", "scratch"),
                   help="Path to scratch dir.")
    p.add_argument("--output-prefix", help="Prefix for output files.")
    p.add_argument("--proportion-unexpressed", 
                   default=0.8, 
                   type=float, 
                   help="Threshold for dropping genes expressed in too few samples.")
    p.add_argument("--variance-filter-level", 
                   default=0.2, 
                   type=float, 
                   help="Threshold for what proportion of the least varying genes to drop.")
    args = p.parse_args()

    # Use input basename as prefix if none was supplied
    prefix = args.output_prefix or os.path.splitext(os.path.basename(args.input_path))[0]

    # Output file path
    normalized_samples_path = os.path.join(args.scratch, "{}log2-normalized.rds".format(prefix))

    expression = prepare_expression(args.input_path, 
                                    args.scratch, 
                                    normalized_samples_path,
                                    prefix=prefix, 
                                    verbose=args.verbose,
                                    proportion_unexpressed=args.proportion_unexpressed, 
                                    variance_filter_level=args.variance_filter_level)
    calculate_correlation(expression, 
                          args.scratch, 
                          prefix=prefix, 
                          verbose=args.verbose)

if __name__ == "__main__":
    main()
