"""
Evaluate Classifier Predictions
Modified from PDX PPTC Machine Learning Analysis
https://github.com/marislab/pdx-classification
Rokita et al. Cell Reports. 2019.
https://doi.org/10.1016/j.celrep.2019.09.071

Gregory Way, 2018
Modified by Krutika Gaonkar for OpenPBTA, 2020

This script evaluates the predictions made by the NF1 and TP53 classifiers in the input PPTC RNAseq data.

## Procedure

1. Load status matrices
  * These files store the mutation status for TP53 and NF1 for the input samples (see 00-tp53-nf1-alterations.R for more information)
2. Align identifiers
  * The identifiers matching the RNAseq data to the status matrix are not aligned.
  * I use an intermediate dictionary to map common identifiers
3. Load predictions (see 01-apply-classifier.py for more details)
4. Evaluate predictions
  * I visualize the distribution of predictions between wild-type and mutant samples for both classifiers

## Output

The output of this notebook are several evaluation figures demonstrating the predictive performance on the input data for the three classifiers.
"""

import os
import random
from decimal import Decimal
from scipy.stats import ttest_ind
import numpy as np
import pandas as pd

from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.metrics import roc_curve, precision_recall_curve

import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] arguments")
parser.add_option(
    "-s", "--statusfile", dest="status_file", help="TP53 and NF1 status file"
)
parser.add_option("-f", "--file", dest="filename", help="scores output file ")
parser.add_option(
    "-c", "--clinical", dest="clinical", help="pbta-histologies.tsv clinical file"
)
parser.add_option(
    "-o",
    "--output_basename",
    dest="outputfile",
    help="output plots basename for TP53 and NF1 ROC curves",
)

(options, args) = parser.parse_args()
status_file = options.status_file
scores_file = options.filename
clinical = options.clinical
outputfilename = options.outputfile


np.random.seed(123)

# read TP53/NF1 alterations
full_status_df = pd.read_table(status_file, low_memory=False)
# read in clinical file
clinical_df = pd.read_table(clinical)
# select only IDs
clinical_df = clinical_df[
    ["Kids_First_Biospecimen_ID", "sample_id", "Kids_First_Participant_ID"]
]

# Obtain a binary status matrix
for idx, val in enumerate(full_status_df.itertuples()):
    # if both SNV/CNV counts are 1 or more than 1 this sample can be considered as a TP53 double hit
    # if SNV or CNV has 1 hit with Li-Fraumeni syndrome cancer predisposition can be considered as double hit
    # if SNV has more than 1 hit sample can be considered as double hit
    if full_status_df.loc[idx, 'SNV_indel_counts'] >= 1 and full_status_df.loc[idx,'CNV_loss_counts'] >= 1:
        full_status_df.loc[idx,'tp53_status'] = 1
    elif full_status_df.loc[idx,'SNV_indel_counts'] == 1 and full_status_df.loc[idx,'cancer_predispositions'] == "Li-Fraumeni syndrome":
        full_status_df.loc[idx,'tp53_status'] = 1
    elif full_status_df.loc[idx,'CNV_loss_counts'] == 1 and full_status_df.loc[idx,'cancer_predispositions'] == "Li-Fraumeni syndrome":
        full_status_df.loc[idx,'tp53_status'] = 1
    elif full_status_df.loc[idx,'SNV_indel_counts'] > 1:
        full_status_df.loc[idx,'tp53_status'] = 1
    else:
        full_status_df.loc[idx,'tp53_status']=0


print("drop tp53_score columns from tp53 annotation file")
full_status_df = full_status_df.drop("tp53_score" , axis = "columns")



# read in scores from 01
file = os.path.join(scores_file)
scores_df = pd.read_table(file)
scores_df = scores_df.rename(str.upper, axis="columns")

scores_df = scores_df.merge(
    full_status_df,
    how="left",
    left_on="SAMPLE_ID",
    right_on="Kids_First_Biospecimen_ID_RNA",
)

print("scores df shape")
print(scores_df.shape)
scores_df.tp53_status.value_counts()

scores_df = scores_df.assign(SAMPLE_ID=scores_df.loc[:, "sample_id"])


gene_status = ["tp53_status"]
# binary counts for tp53 loss status
print("TP53 status")
print(scores_df.tp53_status.value_counts())

print(scores_df.head())

def get_roc_plot(scores_df, gene, outputfilename, color):
    """
    Show roc plot of classifier scores per gene

    Arguments:
    df - the dataframe of scores
    gene - the name of the gene to input
    outputfilename - the name of <filename>_ROC_plot.pdf 
    
    """
    lower_gene = gene.lower()
    scores_df = scores_df.rename(str.lower, axis="columns")
    # Obtain Metrics
    sample_status = scores_df.loc[:, "{}_status".format(lower_gene)]
    sample_score = scores_df.loc[:, "{}_score".format(lower_gene)]
    shuffle_score = scores_df.loc[:, "{}_shuffle".format(lower_gene)]
    print(sample_status.head())
    print(sample_score.head())
    fpr_pbta, tpr_pbta, thresh_pbta = roc_curve(
        sample_status, sample_score, drop_intermediate=False
    )
    precision_pbta, recall_pbta, _ = precision_recall_curve(sample_status, sample_score)
    auroc_pbta = roc_auc_score(sample_status, sample_score)
    aupr_pbta = average_precision_score(sample_status, sample_score)

    # Obtain Shuffled Metrics
    fpr_shuff, tpr_shuff, thresh_shuff = roc_curve(
        sample_status, shuffle_score, drop_intermediate=False
    )
    precision_shuff, recall_shuff, _ = precision_recall_curve(
        sample_status, shuffle_score
    )
    auroc_shuff = roc_auc_score(sample_status, shuffle_score)
    aupr_shuff = average_precision_score(sample_status, shuffle_score)

    roc_df = (
        pd.DataFrame([fpr_pbta, tpr_pbta, thresh_pbta], index=["fpr", "tpr", "threshold"])
        .transpose()
        .assign(gene=gene, shuffled=False, auroc=auroc_pbta)
    )
    # save the dataframe for plotting in R
    roc_df.to_csv(os.path.join("results", outputfilename + "_" + gene + "_roc_threshold_results.tsv"), sep="\t", index=False)
    
    # save shuffled data
    roc_shuff_df = (
        pd.DataFrame([fpr_shuff, tpr_shuff, thresh_shuff], index=["fpr", "tpr", "threshold"])
        .transpose()
        .assign(gene=gene, shuffled=True, auroc=auroc_shuff)
    )
    
    roc_shuff_df.to_csv(os.path.join("results", outputfilename + "_" + gene + "_roc_threshold_results_shuffled.tsv"), sep="\t", index=False)

    plt.subplots(figsize=(5, 5))
    plt.axis("equal")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.plot(
        fpr_pbta,
        tpr_pbta,
        label="{} (AUROC = {})".format(gene, round(auroc_pbta, 2)),
        linestyle="solid",
        color=color,
    )

    # Shuffled Data
    plt.plot(
        fpr_shuff,
        tpr_shuff,
        label="{} Shuffle (AUROC = {})".format(gene, round(auroc_shuff, 2)),
        linestyle="dotted",
        color=color,
    )

    plt.xlabel("False Positive Rate", fontsize=12)
    plt.ylabel("True Positive Rate", fontsize=12)
    plt.tick_params(labelsize=10)

    lgd = plt.legend(bbox_to_anchor=(0.3, 0.15), loc=2, borderaxespad=0.0, fontsize=10)
    plt.savefig(os.path.join("plots", outputfilename + "_" + gene + "_roc.png"))

get_roc_plot(scores_df, gene="TP53", outputfilename=outputfilename, color="#7570b3")
