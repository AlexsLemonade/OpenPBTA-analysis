
"""
Filter polyA and ribodeplete (stranded) samples to only tumor samples that are MEND QC pass.
Then, create correlation matrices for these samples using gene expression data.
The correlations will be calculated using pairwise Spearman correlation of the RNA-Seq gene expression profiles.
The gene expression profile data will be filtered using a developed method described in Vaske et al.

https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/229
"""

import os
import sys
import argparse
import math
import tarfile
import numpy as np
import pandas as pd
import scipy.stats
import sklearn.metrics.pairwise as sklp
import utils

def get_qc_status(mend_results_path, mend_manifest_path):
    """Takes: mend results tgz containing bam_umend_qc.tsv files, and
    manifest containing mapping from filename to sample id.
    returns dictionary from sample id to mend status: PASS or FAIL"""
    mend_tgz = tarfile.open(mend_results_path)
    qc_files = (i for i in mend_tgz.getmembers() if i.name.endswith("bam_umend_qc.tsv"))
    # Dictionary of filename (UUID.bam_umend_qc.tsv) to PASS or FAIL string
    filename_map = { i.name : mend_tgz.extractfile(i).readlines()[1][-5:-1].decode("utf-8") for i in qc_files }
    # Find sample ID entries in manifest and map them via filename.
    manifest = utils.read_tsv(mend_manifest_path)
    return { manifest.loc[k]["Kids.First.Biospecimen.ID"] : v for k, v in filename_map.items() }

def filter_samples(input_matrix,
                   mend_results_path,
                   mend_manifest_path,
                   clinical_path,
                   nofilter=False,
                   verbose=False):
    """Filters samples by QC and tumor status.
    Takes path to samples expression file; tgz of MEND QC results;
    manifest of MEND QC files to sample ids; and clinical data file.
    Returns data frame of sample expression filtered to contain only
    MEND QC pass tumor samples. Filters are skippable via nofilter if,
    for example, clinical and qc files aren't available."""
    print_v = print if verbose else lambda *a, **k: None

    print_v("Loading sample matrix {}".format(input_matrix))
    raw_samples = utils.read_rds(input_matrix)

    if nofilter:
        print_v("Skipped the QC and tumor status filters; using all samples.")
        return raw_samples

    # Filter to only samples with a MEND QC status of PASS
    qc_status = get_qc_status(mend_results_path, mend_manifest_path)
    qc_pass_ids= [k for k, v in qc_status.items() if v == "PASS"]
    qc_pass_samples = raw_samples[set(raw_samples.columns).intersection(qc_pass_ids)]
    print_v("QC filter: dropped {} samples that were not MEND QC PASS.".format(
            len(raw_samples.columns) - len(qc_pass_samples.columns)))

    # Further filter to only tumor RNA-Seq samples
    clinical = utils.read_tsv(clinical_path)
    tumor_ids = clinical[(clinical["sample_type"] == "Tumor")
                         & (clinical["composition"] == "Solid Tissue")
                         & (clinical["experimental_strategy"] == "RNA-Seq")
                        ].index
    qc_pass_tumor_samples = qc_pass_samples[set(qc_pass_samples.columns).intersection(tumor_ids)]
    print_v("Tumor filter: dropped {} samples that were not solid tumor RNA-Seq.".format(
            len(qc_pass_samples.columns) - len(qc_pass_tumor_samples.columns)))
    return qc_pass_tumor_samples

def prepare_expression(raw_samples,
                       normalized_samples_path,
                       gene_filters_path,
                       proportion_unexpressed=0.8,
                       variance_filter_level=0.2,
                       verbose=False):
    """Converts sample matrix to log2(TPM+1) and applies expression and variance filters.
    returns dataframe of converted and filtered expression; generates output feather files
    stored at normalized_samples_path, gene_filters_path.
    raw_samples = TPM expression dataframe.
    Normalized samples file is a log2(TPM+1) expression matrix retaining all genes.
    Gene filters file lists all genes and their filter status - either R = Retained,
    E = dropped by the Expression filter, or V = dropped by the variance filter."""

    print_v = print if verbose else lambda *a, **k: None

    # Convert to log2(tpm+1)
    samples = np.log2(raw_samples+1)
    print_v("Writing normalized samples to {}".format(normalized_samples_path))
    utils.write_feather(samples, normalized_samples_path)

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

    # Save filter status by first by marking every gene as failed Expression filter;
    # marking genes that passed as failed Variance filter, and again marking genes that passed
    # that filter as Retained.
    gene_filter_status = pd.DataFrame("E", index=samples.index, columns=["status"])
    gene_filter_status["status"][expression_filtered_df.index] = "V"
    gene_filter_status["status"][expr_var_filtered_df.index] = "R"

    print_v("Writing gene filter status to {}".format(gene_filters_path))
    utils.write_feather(gene_filter_status, gene_filters_path)

    return expr_var_filtered_df

def calculate_correlation(expr_var_filtered_df, all_by_all_path, verbose=False):
    """Takes an expression matrix dataframe and calculates spearman correlations.
    (spearman is a rank transformed pearson ('correlation') metric)
    output file: scratch/{prefix}all_by_all_correlations.feather"""
    print_v = print if verbose else lambda *a, **k: None

    # Rank transform - values are now that rank of that gene within each sample
    print_v("Rank transforming")
    filtered_rank_transformed_df = np.apply_along_axis(scipy.stats.rankdata, 0, expr_var_filtered_df)

    # Run pearson metric on ranked genes. For pairwise_distances we need samples=rows so transpose.
    print_v("Running pairwise distances")
    x_corr = 1 - sklp.pairwise_distances(X=filtered_rank_transformed_df.transpose(), metric="correlation")

    # Reapply the axis labels
    labels=expr_var_filtered_df.columns
    all_by_all_df = pd.DataFrame(x_corr, index=labels, columns=labels)

    print_v("Writing to file {}".format(all_by_all_path))
    utils.write_feather(all_by_all_df, all_by_all_path)
    return all_by_all_df

def main():
    """Creates correlation matrices for dataset."""
    # This script should always run as if it were being called from
    # the directory it lives in.
    os.chdir(sys.path[0])

    p = argparse.ArgumentParser()
    p.add_argument("input_path", metavar="input-path", help="Path to input Rds file.")
    p.add_argument("--clinical-path", help="Path to clinical data matrix tsv file.")
    p.add_argument("--qc-manifest-path", help="Path to MEND QC manifest tsv file.")
    p.add_argument("--qc-results-path", help="Path to tgz file of MEND QC results for each sample.")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--nofilter", action="store_true", help="Skip the sample QC and tumor status filters")
    p.add_argument("--scratch", 
                   default=os.path.join("..", "..", "scratch"),
                   help="Path to scratch dir.")
    p.add_argument("--prefix", help="Prefix for output files.")
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
    prefix = args.prefix or os.path.splitext(os.path.basename(args.input_path))[0]

    # Output file paths
    normalized_samples_path = os.path.join(args.scratch, "{}log2-normalized.feather".format(prefix))
    gene_filters_path = os.path.join(args.scratch, "{}gene-filters.feather".format(prefix))
    all_by_all_path = os.path.join(args.scratch, "{}all_by_all_correlations.feather".format(prefix))


    # Load the expression file and filter samples to QC-pass tumors only
    filtered_expression = filter_samples(args.input_path,
                                         args.qc_results_path,
                                         args.qc_manifest_path,
                                         args.clinical_path,
                                         nofilter=args.nofilter,
                                         verbose=args.verbose)

    normalized_expression = prepare_expression(filtered_expression,
                                               normalized_samples_path,
                                               gene_filters_path,
                                               proportion_unexpressed=args.proportion_unexpressed,
                                               variance_filter_level=args.variance_filter_level,
                                               verbose=args.verbose)
    calculate_correlation(normalized_expression,
                          all_by_all_path,
                          verbose=args.verbose)

if __name__ == "__main__":
    main()
