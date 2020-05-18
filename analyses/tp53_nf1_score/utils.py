"""
Import helper methods
Gregory Way 2018

Usage:

    from utils import apply_classifier
    from utils import shuffle_columns
"""


def apply_classifier(coef_df, rnaseq_df):
    """
    Apply the logistic regression classifier coefficients to the input gene
    expression matrix. The function will perform the gene subsetting.

    Arguments:
    coef_df - a pandas dataframe with feature and feature weight (importance)
    rnaseq_df - a scaled pandas dataframe (sample by gene)

    Output:
    A tuple of classifier scores for the input data, the genes in common
    between the two inputs, and the missing genes
    """

    import numpy as np

    rnaseq_genes = set(rnaseq_df.columns)
    classifier_genes = set(coef_df["feature"])

    # Determine the extent of coefficient overlap
    common_coef = list(classifier_genes & rnaseq_genes)
    common_coef = coef_df.query("feature in @common_coef")

    # Which Genes are Missing?
    missing_genes = list(classifier_genes.difference(rnaseq_genes))
    missing_genes = coef_df.query("feature in @missing_genes")

    # Get the weights ready for applying the classifier
    apply_weights = common_coef.loc[:, ["feature", "weight"]]
    apply_weights = apply_weights.set_index("feature")

    # Make sure the indeces match between RNAseq and coefficients
    # This will also subset the RNAseq matrix
    rnaseq_df = rnaseq_df.reindex(apply_weights.index, axis="columns")

    # Apply a logit transform [y = 1/(1+e^(-wX))] to output scores
    scores = apply_weights.T.dot(rnaseq_df.T)
    scores = 1 / (1 + np.exp(-1 * scores))

    return (scores, common_coef, missing_genes)


def shuffle_columns(gene):
    """
    To be used in an `apply` pandas func to shuffle columns around a datafame
    Import only
    """
    import numpy as np

    return np.random.permutation(gene.tolist())
