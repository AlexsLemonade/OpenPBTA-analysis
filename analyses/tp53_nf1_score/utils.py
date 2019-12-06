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
    classifier_genes = set(coef_df['feature'])

    # Determine the extent of coefficient overlap
    common_coef = list(classifier_genes & rnaseq_genes)
    common_coef = coef_df.query('feature in @common_coef')

    # Which Genes are Missing?
    missing_genes = list(classifier_genes.difference(rnaseq_genes))
    missing_genes = coef_df.query('feature in @missing_genes')

    # Get the weights ready for applying the classifier
    apply_weights = common_coef.loc[:, ['feature', 'weight']]
    apply_weights = apply_weights.set_index('feature')

    # Make sure the indeces match between RNAseq and coefficients
    # This will also subset the RNAseq matrix
    rnaseq_df = rnaseq_df.reindex(apply_weights.index, axis='columns')

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


def perform_ttest(df, gene):
    """
    Perform an independent ttest for the status predictions for each gene class

    Arguments:
    df - the dataframe of scores
    gene - the name of the gene to subset the scores dataframe
    """
    from scipy.stats import ttest_ind

    lower_gene = gene.lower()

    if lower_gene not in ['ras', 'nf1', 'tp53']:
        raise ValueError('Enter either "ras", "nf1", or "tp53" as the gene')

    status = '{}_status'.format(lower_gene)
    score = '{}_score'.format(lower_gene)

    one_status = '{} == 1'.format(status)
    zero_status = '{} == 0'.format(status)

    t_results = ttest_ind(
        a=df.query(one_status).loc[:, score],
        b=df.query(zero_status).loc[:, score],
        equal_var=False
        )

    return t_results


def get_mutant_boxplot(df, gene, t_test_results=None, histology=False,
                       hist_color_dict=None):
    """
    Show boxplot of classifier scores stratified by binary gene mutation status

    Arguments:
    df - the dataframe of scores
    gene - the name of the gene to input
    t_test_results - the output of `perform_ttest` for the specific gene
    histology - boolean if the boxplot should stratify by histology
    hist_color_dict - a dictionary storing prespecified hex colors by histology
    """
    import os
    from decimal import Decimal
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Set plotting constants
    x1, x2 = 0, 1
    y1, h = 0.98, 0.03

    lower_gene = gene.lower()

    if lower_gene not in ['ras', 'nf1', 'tp53']:
        raise ValueError('Enter either "ras", "nf1", or "tp53" as the gene')

    x = '{}_status'.format(lower_gene)
    y = '{}_score'.format(lower_gene)
    xticklabels = ['{} Wild-Type'.format(gene), '{} Mutant'.format(gene)]

    if histology:
        output_file = os.path.join('figures',
                                   'histology_{}_predictions.pdf'.format(gene))
        plt.rcParams['figure.figsize'] = (9, 4)
        ax = sns.boxplot(x=x,
                         y=y,
                         data=df,
                         hue='Histology.Detailed',
                         color='white',
                         fliersize=0)

        ax = sns.stripplot(x=x,
                           y=y,
                           data=df,
                           hue='Histology.Detailed',
                           palette=hist_color_dict,
                           dodge=True,
                           edgecolor='black',
                           jitter=0.25,
                           size=4,
                           alpha=0.65)

        handles, labels = ax.get_legend_handles_labels()

        lgd = plt.legend(handles[25:50], labels[25:50],
                         ncol=2,
                         bbox_to_anchor=(1.03, 1),
                         loc=2,
                         borderaxespad=0.,
                         fontsize=8)
        lgd.set_title("Histology")

    else:
        output_file = os.path.join('figures',
                                   '{}_predictions.pdf'.format(gene))
        plt.rcParams['figure.figsize'] = (3.5, 4)
        ax = sns.boxplot(x=x,
                         y=y,
                         data=df,
                         palette='Greys',
                         fliersize=0)

        ax = sns.stripplot(x=x,
                           y=y,
                           data=df,
                           dodge=True,
                           edgecolor='black',
                           jitter=0.25,
                           size=4,
                           alpha=0.65)

    ax.set_ylabel('Classifier Score', fontsize=12)
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_xticklabels(xticklabels)

    # Add Ras T-Test Results
    if not histology:
        plt.plot([x1, x1, x2, x2], [y1, y1+h, y1+h, y1], lw=1.2, c='black')
        plt.text(.5, y1+h, "{:.2E}".format(Decimal(t_test_results.pvalue)),
                 ha='center', va='bottom', color="black")

    plt.axhline(linewidth=2, y=0.5, color='black', linestyle='dashed')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


def vis_classifier_scores(df, gene, rcparam=(6, 4), variant_plot=False):
    """
    Build boxplot of classifier scores stratified by gene confidence scores

    Arguments:
    df - the dataframe of scores
    gene - the name of the gene to input
    variant_plot - boolean if the x axis should be variant classifications
    """
    import os
    import matplotlib.pyplot as plt
    import seaborn as sns

    lower_gene = gene.lower()

    if lower_gene not in ['ras', 'nf1', 'tp53']:
        raise ValueError('Enter either "Ras", "NF1", or "TP53" as the gene')

    y = '{}_score'.format(lower_gene)
    ylabel = '{} Classifier Score'.format(gene)
    output_file = os.path.join('figures', '{}_scores.pdf'.format(gene))

    plt.rcParams['figure.figsize'] = rcparam
    ax = sns.boxplot(x='Hugo_Symbol',
                     y=y,
                     data=df,
                     color='white',
                     fliersize=0)

    ax = sns.stripplot(x='Hugo_Symbol',
                       y=y,
                       data=df,
                       dodge=True,
                       edgecolor='black',
                       jitter=0.25,
                       size=4,
                       alpha=0.65)

    ax.set_ylim([0, 1])
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel('Genes', fontsize=12)
    handles, labels = ax.get_legend_handles_labels()

    plt.legend(handles[4:8], labels[4:8],
               ncol=1,
               bbox_to_anchor=(1.03, 1),
               loc=2,
               borderaxespad=0.,
               fontsize=10)

    plt.axhline(linewidth=2, y=0.5, color='black', linestyle='dashed')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


def extract_outliers(df, gene):
    """
    Extract positive and negative outliers in score dataframes

    Arguments:
    df - the gene specific scores dataframe
    gene - the gene of interest

    Output:
    A summary dataframe of samples predicted incorrectly
    """
    import pandas as pd

    lower_gene = gene.lower()

    if lower_gene not in ['ras', 'nf1', 'tp53']:
        raise ValueError('Enter either "Ras", "NF1", or "TP53" as the gene')

    score = '{}_score'.format(lower_gene)

    # Obtain false positives and false negatives
    false_negatives = (
        df
        .query('Hugo_Symbol != "wild-type"')
        .sort_values(by=score, ascending=False)
        .iloc[:, 0:3]
    )
    false_negatives = false_negatives[false_negatives[score] < 0.5]

    false_positives = (
        df
        .query('Hugo_Symbol == "wild-type"')
        .sort_values(by=score, ascending=False)
        .iloc[:, 0:3]
    )
    false_positives = false_positives[false_positives[score] >= 0.5]

    col_names = ['sample_id', 'classifier_score', 'histology']
    false_negatives.columns = col_names
    false_positives.columns = col_names

    false_negatives = false_negatives.assign(hugo_symbol=gene,
                                             true_status='altered')
    false_positives = false_positives.assign(hugo_symbol=gene,
                                             true_status='wild-type')

    # Combine incorrect predictions and output
    outliers_df = pd.concat([false_negatives, false_positives], axis='rows')
    return outliers_df.sort_values(by='classifier_score', ascending=False)
