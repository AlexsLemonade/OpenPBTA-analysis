"""
Generates gene outlier thresholds for normalized samples and uses them to
calculate the outliers.

Outlier thresholds are calculated according to a modified Tukey method, detailed
below.

Then, outliers are calculated based on these thresholds. In addition, for each sample,
the 5% of genes with the highest expression are marked, as well as the genes which
are excluded by the expression and variance filters from step 1.

The final result is a matrix .rds file of genes by samples where, for each sample,
each gene is labeled with a 4-character string formatted as follows:

"(P|F) (U|D|N) (T|L) (R|E|V)"

0. QC status of this sample - do at least 5% of genes in this sample have
expression greater than zero?
P: Pass, F: Fail

1. Outlier status:
U: Up outlier, D: Down outlier, N: Not an outlier

2. Is this gene's expression in the top 5 percent of expression for this sample?
T: in Top 5%, L: in Lower 95%

3. Is this gene dropped from the analysis due to expression or variance filters
(as implemented in step 1)?
R: Retained, E: dropped by Expression filter, V: dropped by Variance filter
(TODO: E/V distinction is not yet implemented -- all dropped genes are labeled E.)

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
                        iqr_multiplier=1.5,
                        verbose=False):
    """Generate high and low outlier thresholds for each gene
    using the modified tukey method from a matrix of
    log2(TMP+1) normalized sample expression"""
    print_v = print if verbose else lambda *a, **k: None
    print_v("Calculating outlier thresholds.")
    return normalized_samples.apply(modified_tukey_thresholds, args=[iqr_multiplier], axis=1)

def update_at_index(idx, new_value):
    """Helper function to modify list in-place, setting value at idx to new_value.
    for use in pandas apply and applymap, eg x.apply(update_at_index(0, "A"))."""
    def updater(mylist):
        mylist[idx] = new_value
    return updater

def single_sample_outliers(sample_expression, expression_thresholds):
    """Given the expression of a single sample as series, and matrix
    of thresholds, calculate the up and down outliers and genes in the top 5%
    of sample's expression. Return a matrix of genes vs samples, with each gene
    marked with, in order, QC status(Pass/Fail), outlier status (Up/Down/Not outlier),
    and top5% status(Top5/Lower95); additionally mark all genes as Retained by filter, 
    to be modified by parent."""

    down_outliers = sample_expression[sample_expression < expression_thresholds["low"]].index
    up_outliers = sample_expression[sample_expression > expression_thresholds["high"]].index

    # Mark genes with expression in the top 5% within that sample (inclusive)
    top5_value = sample_expression.quantile(0.95)
    qc_fail = (top5_value <=0) # If 0 is in the top5%, this sample is really underexpressed
    top5_genes = sample_expression[sample_expression >= top5_value].keys()

    # Generate result code. Start with default Pass/Not outlier/Lower95/Retained.
    result = pd.Series("PNLR", index = sample_expression.index)
    result = result.apply(list)# Split into a list to allow modification by index
    if qc_fail:
        result.apply(update_at_index(0, "F")) # Mark every gene as QC fail for this sample
    result[up_outliers].apply(update_at_index(1, "U"))
    result[down_outliers].apply(update_at_index(1, "D"))
    result[top5_genes].apply(update_at_index(2, "T"))
    return result


def calculate_all_outliers(samples,
                           thresholds,
                           filtered_genelist_path,
                           outliers_path,
                           verbose=False):
    """For each sample calculate up and down outlier genes based on the thresholds.
    Then, apply expression & variance filters to mark 'noise' genes, and save as result file."""
    print_v = print if verbose else lambda *a, **k: None
    print_v("Calculating outliers for all samples")

    outliers = samples.apply(single_sample_outliers, args=[thresholds], axis=0)

    # Apply the expression & variance filters that we calculated in step 1 to mark filtered
    # genes as dropped due to Expression or Variance; by default they're already marked Retained.
    filtered_genelist = utils.read_rds(filtered_genelist_path)
    dropped_genes = pd.Index(set(outliers.index) - set(filtered_genelist["gene_id"]))
    outliers.loc[dropped_genes].applymap(update_at_index(3, "E"))
    outliers = outliers.applymap("".join) # Concatenate all values to string

    print_v("Writing outlier results to {}".format(outliers_path))
    utils.write_tsv(outliers, outliers_path)
    return outliers

def main():
    """Generates expression thresholds and uses them to calculate
    high and low outliers for each sample in the dataset.
    input file: in format log_2(TPM+1), samples in columns, genes in rows
    output file: genes are rows, samples are columns,
    values are a fixed width string of keys representing
    QC, outlier, top5, and dropped-gene status."""
    # This script should always run as if it were being called from
    # the directory it lives in.
    os.chdir(sys.path[0])

    p = argparse.ArgumentParser()
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--scratch",
                   default=os.path.join("..", "..", "scratch"),
                   help="Path to scratch dir.")
    p.add_argument("--results",
                   default="results",
                   help="Path to results dir for generated outliers file.")
    p.add_argument("--prefix", help="Prefix for input and output files.")
    p.add_argument("--iqr-multiplier",
                   default=1.5,
                   type=float,
                   help="Interquartile range multiplier for Tukey outlier calculation.")
    args = p.parse_args()

    print_v = print if args.verbose else lambda *a, **k: None

    # Input paths
    normalized_samples_path = os.path.join(args.scratch, "{}log2-normalized.rds".format(args.prefix))
    filtered_genelist_path = os.path.join(args.scratch,
                                          "{}filtered_genes_to_keep.rds".format(args.prefix))
    # Output paths
    outliers_path = os.path.join(args.results, "{}gene_expression_outliers.tsv.gz".format(args.prefix))

    os.makedirs(args.results, exist_ok=True) # make results dir if not already present

    # Load normalized samples from Step 01
    normalized_samples = utils.read_rds(normalized_samples_path)

    # Generate thresholds and outliers
    thresholds = generate_all_thresholds(normalized_samples,
                                     iqr_multiplier=args.iqr_multiplier,
                                     verbose=args.verbose)
    calculate_all_outliers(normalized_samples,
                           thresholds,
                           filtered_genelist_path,
                           outliers_path,
                           verbose=args.verbose)
    print_v("Done generating thresholds and outliers.")

if __name__ == "__main__":
    main()
