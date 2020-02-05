"""
Generates gene outlier thresholds for normalized samples and uses them to
calculate the outliers.

Outlier thresholds are calculated according to a modified Tukey method, detailed
below.

Then, outliers are calculated based on these thresholds. In addition, for each sample,
the 5% of genes with the highest expression are marked, as well as the genes which
are excluded by the expression and variance filters from step 1.

The final result is a matrix .rds file of genes by samples where, for each sample,
each gene may be labeled with zero or more of the following keys (concatenated
together into a string):

U: Up outlier
D: Down outlier
T: in Top five percent of genes for this sample
F: qc Fail (when less than 5% of genes in the sample have nonzero expression)
P: gene is droPped by expression or variance filter


Modified Tukey Method
--------------------

The current sample being considered for outlier expression is dropped from the
dataset, and then the Tukey outlier thresholds of
first_quartile - (1.5 * IQR) and third_quartile + (1.5 * IQR) are calculated
for each gene. This allows outlier results to be concordant whether a focus
sample was originally present in the dataset or whether it is a new sample.

Running this calculation for each sample in the dataset is computationally
expensive; but there is a modification we can make which allows us to generate
a single set of thresholds to be used on all samples.

Consider a focus sample s that for gene g has expression Esg. Removing this
sample from the dataset, we generate a personalized threshold Tsg for that gene.
What we'd like to generate is a generalized threshold Tg that is valid for
all samples. Tg doesn't need to be equal to Tsg, but it must always return
the same outlier status, ie, one of the following must be true:

Esg < Tsg < Tg
Esg < Tg < Tsg
Tsg < Tg < Esg
Tg < Tsg < Esg

The key observations are that the values of the quartiles and IQR depend only
on the number of observations and on the values of those observations directly
at the quartile boundaries; and that any outlier values will be far out into quartiles
1 or 4 and thus the quartile + IQR values can't possibly depend directly on their value.
Thus, if for gene g we drop any single expression value that is an up outlier, the
subsequently calculated threshold will be correct for all samples for which g is
an up outliers or almost-up-outlier. The threshold will be slightly incorrect for
values that are at or below the 3rd quantile, but those are nowhere near being
up outliers so that doesn't matter. The same principle holds for down outliers.

Therefore, we modify the Tukey method by, for each gene, dropping the single highest
expression value (for simplicity) and calculating the high threshold; and similarly
dropping the single lowest expression value and calculating the low threshold.

https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/229
"""

import os
import math
import sys
import argparse
import pandas as pd
import numpy as np
import utils

def normalize_expression(input_path, scratch_dir, prefix="", verbose=False):
    """Convert samples to log2(TPM) + 1 and save in scratch.
    Note: This is duplicate work from step 01 because I forgot in that step
    to write results to scratch file. Future update will merge these. (TODO)"""
    print_v = print if verbose else lambda *a, **k: None
    print_v("Loading samples {}".format(input_path))
    raw_samples = utils.read_rds(input_path)
    samples = np.log2(raw_samples+1)
    output_path = os.path.join(scratch_dir, "{}log2-normalized.rds".format(prefix))
    print_v("Writing normalized samples to {}".format(output_path))
    utils.write_rds(samples, output_path)
    return samples

def tukey_thresholds(expression_values, iqr_multiplier):
    """Generate a pair of outlier thresholds (high, low) from
    expression_values (series) and an iqr multiplier using standard Tukey method."""
    quartile1 = np.percentile(expression_values, 25)
    quartile3 = np.percentile(expression_values, 75)
    iqr = quartile3 - quartile1
    outlier_range = iqr_multiplier * iqr
    return (quartile3 + outlier_range, quartile1 - outlier_range)

def modified_tukey_thresholds(expression_values, iqr_multiplier):
    """Use a modified Tukey method to generate outlier thresholds (high, low)
    for a single gene. Drops the highest value when calculating the high threshold,
    and the lowest value for low threshold, to match thresholds for outlier
    status when a sample is not present in the original dataset."""
    drop_highest = expression_values.drop(expression_values.idxmax())
    high = tukey_thresholds(drop_highest, iqr_multiplier)[0]
    drop_lowest = expression_values.drop(expression_values.idxmin())
    low = tukey_thresholds(drop_lowest, iqr_multiplier)[1]
    return( pd.Series([high, low], index=["high", "low"]))

def generate_all_thresholds(normalized_samples,
                        scratch_dir,
                        iqr_multiplier=1.5,
                        prefix="",
                        verbose=False):
    """Generate high and low outlier thresholds for each gene
    using the modified tukey method from a matrix of
    log2(TMP+1) normalized sample expression"""
    print_v = print if verbose else lambda *a, **k: None

    print_v("Calculating outlier thresholds.")
    thresholds = normalized_samples.apply(modified_tukey_thresholds, args=[iqr_multiplier], axis=1)
    thresholds_path = os.path.join(scratch_dir, "{}threshold-expression-values.rds".format(prefix))
    print_v("Writing expression outlier thresholds to {}".format(thresholds_path))
    utils.write_rds(thresholds, thresholds_path)
    return thresholds

def single_sample_outliers(sample_expression, expression_thresholds):
    """Given the expression of a single sample as series, and matrix
    of thresholds, calculate the up and down outliers and genes in the top 5%
    of sample's expression. Return a matrix of genes vs samples, with each gene
    marked as U or D for outlier, T if it's in the top 5%, and F for qc Fail."""

    down_outliers = sample_expression[sample_expression < expression_thresholds["low"]].index
    up_outliers = sample_expression[sample_expression > expression_thresholds["high"]].index

    # Mark genes with expression in the top 5% within that sample (inclusive)
    top5_value = sample_expression.quantile(0.95)
    qc_fail = (top5_value <=0) # If 0 is in the top5%, this sample is really underexpressed
    top5_genes = sample_expression[sample_expression >= top5_value].keys()

    # Possible values include (D)own outlier, (U)p outlier, (T)op5%, and qc (F)ail
    # dro(P)ped genes are added in calculate_outliers
    result = pd.Series("", index = sample_expression.index)
    result[top5_genes] += "T"
    result[down_outliers] += "D"
    result[up_outliers] += "U"
    if qc_fail:
        result += "F" # Mark every gene as QC fail for this sample
    return result

def calculate_all_outliers(samples,
                           thresholds,
                           scratch_dir,
                           results_dir,
                           prefix="",
                           verbose=False):
    """For each sample calculate up and down outlier genes based on the thresholds.
    Then, apply expression & variance filters to mark 'noise' genes, and save as result file."""
    print_v = print if verbose else lambda *a, **k: None
    print_v("Calculating outliers for all samples")

    outliers = samples.apply(single_sample_outliers, args=[thresholds], axis=0)

    # Apply the expression & variance filters that we calculated in step 1 to mark filtered
    # genes as dro(P)ped.
    filtered_genelist_filepath = os.path.join(scratch_dir,
                                              "{}filtered_genes_to_keep.rds".format(prefix))
    filtered_genelist = utils.read_rds(filtered_genelist_filepath)
    dropped_genes = pd.Index(set(outliers.index) - set(filtered_genelist["gene_id"]))
    outliers.loc[dropped_genes] = outliers.loc[dropped_genes].applymap(lambda x: x + "P")

    outliers_path = os.path.join(results_dir, "{}gene_expression_outliers.rds".format(prefix))
    print_v("Writing outlier results to {}".format(outliers_path))
    os.makedirs(results_dir, exist_ok=True) # make results dir if not already present
    utils.write_rds(outliers, outliers_path)
    return outliers

def main():
    """Generates expression thresholds and uses them to calculate
    high and low outliers for each sample in the dataset.
    input file: in format TPM, samples in columns, genes in rows
    output file: (TODO finalize): genes are rows, samples are columns,
    values are a comma-separated list of outlier and top5 status, with
    dropped genes removed."""
    # This script should always run as if it were being called from
    # the directory it lives in.
    os.chdir(sys.path[0])

    p = argparse.ArgumentParser()
    p.add_argument("input_path", metavar="input-path",
                   help="Path to input Rds file. Genes in rows, samples in columns, TPM format.")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--scratch",
                   default=os.path.join("..", "..", "scratch"),
                   help="Path to scratch dir.")
    p.add_argument("--results",
                   default="results",
                   help="Path to results dir for generated outliers file.")
    p.add_argument("--output-prefix", help="Prefix for output files.")
    p.add_argument("--iqr-multiplier",
                   default=1.5,
                   type=float,
                   help="Interquartile range multiplier for Tukey outlier calculation.")
    args = p.parse_args()

    print_v = print if args.verbose else lambda *a, **k: None

    # Use input basename as prefix if none was supplied
    prefix = args.output_prefix or os.path.splitext(os.path.basename(args.input_path))[0]

    # Load normalized samples if they exist; otherwise, generate.
    normalized_samples_path = os.path.join(args.scratch, "{}log2-normalized.rds".format(prefix))
    try:
        normalized_samples = utils.read_rds(normalized_samples_path)
        print_v("Loaded existing normalized samples from file.")
    except utils.pyreadr.custom_errors.LibrdataError:
        normalized_samples = normalize_expression(args.input_path,
                                                  args.scratch,
                                                  prefix=prefix,
                                                  verbose=args.verbose)
    # Generate thresholds and outliers
    thresholds = generate_all_thresholds(normalized_samples,
                                     args.scratch,
                                     iqr_multiplier=args.iqr_multiplier,
                                     prefix=prefix,
                                     verbose=args.verbose)
    calculate_all_outliers(normalized_samples,
                           thresholds,
                           args.scratch,
                           args.results,
                           prefix=prefix,
                           verbose=args.verbose)
    print_v("Done generating thresholds and outliers.")

if __name__ == "__main__":
    main()
