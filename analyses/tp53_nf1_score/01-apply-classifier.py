"""
Apply Machine Learning Classifiers to OpenPBTA data
Modified from PDX PPTC Machine Learning Analysis
https://github.com/marislab/pdx-classification
Rokita et al. Cell Reports. 2019. 
https://doi.org/10.1016/j.celrep.2019.09.071

Gregory Way, 2018
Modified by Krutika Gaonkar for OpenPBTA, 2019
 
In the following notebook, three distinct classifiers to OpenPedCan data (TPM).
The first classifier detects Ras activation. For more details about the algorithm and results, refer to [Way et al. 2018](https://doi.org/10.1016/j.celrep.2018.03.046 "Machine Learning Detects Pan-cancer Ras Pathway Activation in The Cancer Genome Atlas"). I also include _TP53_ inactivation predictions. This classifier was previously applied in [Knijnenburg et al. 2018](https://doi.org/10.1016/j.celrep.2018.03.076 "Genomic and Molecular Landscape of DNA Damage Repair Deficiency across The Cancer Genome Atlas"). The third is a classifier that predicts _NF1_ inactivation. We previously applied this in [Way et al. 2017](https://doi.org/10.1186/s12864-017-3519-7 "A machine learning classifier trained on cancer transcriptomes detects NF1 inactivation signal in glioblastoma").
 
To apply other classifiers (targeting other genes of interest) refer to https://github.com/greenelab/pancancer.

## Procedure
 
1. Load RNAseq matrix (`.RDS`)
  * The matrix is in `gene symbol` x `sample` format 
2. Process matrix
  * Take the z-score by gene
3. Load classifier coefficients
  * For both Ras, TP53, and NF1 classifiers
4. Apply each classifier to the input data
  * This also requires additional processing steps to the input data (subsetting to the respective classifier genes)
  * Also note that not all genes are present in the input RNAseq genes.
5. Shuffle the gene expression data (by gene) and apply each classifier to random data
 
### Important Caveat
 
Because not all of the classifier genes are found in the input dataset, the classifier probability is not calibrated correctly. The scores should be interpreted as continuous values representing relative gene alterations, and not as a pure probability estimate.

## Output
 
The output of this notebook are the predicted scores for both classifiers across all samples for real data and shuffled data. This is in the form of a single text file with three columns (`sample_id`, `ras_score`, `tp53_score`, `nf1_score`,  `ras_shuffle`, `tp53_shuffle`, `nf1_shuffle`).
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from sklearn.preprocessing import StandardScaler
from rpy2.robjects import pandas2ri
from utils import apply_classifier, shuffle_columns
from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] arguments")
parser.add_option(
    "-f", "--expfile", dest="expfile", help="rds file sample X genes expression"
)
parser.add_option(
    "-t", "--histology", dest="histology", help="histology file for subsetting expression file"
)
parser.add_option(
     "-c", "--cohorts", type ='string', dest="cohorts", help="list of cohorts of interest"
 )

(options, args) = parser.parse_args()
exprs_file = options.expfile
histology = options.histology
cohort_list = options.cohorts.split(",")

np.random.seed(123)
pandas2ri.activate()
readRDS = robjects.r["readRDS"]
rownamesRDS = robjects.r["rownames"]

# Load gene expression data in rds format
name = os.path.basename(exprs_file)
name = re.sub("\.rds$", "", name)

exprs_rds = readRDS(exprs_file)
exprs_total = pandas2ri.ri2py(exprs_rds)
exprs_total.index = rownamesRDS(exprs_rds)


# Read in histologies file
histology = pd.read_csv(histology, sep="\t")
# filter histology file based on cohort
histology = histology[histology['cohort'].isin(cohort_list)]
histology = histology[histology.cohort != "TCGA"]
histology = histology[histology.cohort != "GTEx"]
cohort_list = list(histology['cohort'].unique())

exprs_scaled_df = pd.DataFrame()
exprs_shuffled_df = pd.DataFrame()

for cohort in cohort_list:
  # Filter the expression_df to tumor, RNA-Seq and cohort
  histology_cohort = histology[(histology["experimental_strategy"] == "RNA-Seq") & (histology["sample_type"] == "Tumor") & (histology["cohort"] == cohort)]
  RNA_libraries = list(set(histology_cohort['RNA_library']))
  
  # Remove NAN from RNA library list
  RNA_libraries = [RNA_libraries for RNA_libraries in RNA_libraries if str(RNA_libraries) != 'nan']

  for rna_library in RNA_libraries:
    tumor_RNA_cohort = list(histology_cohort[histology_cohort["RNA_library"] == rna_library]['Kids_First_Biospecimen_ID'])
    exprs_df = exprs_total[tumor_RNA_cohort]

    # Transpose the expression dataframe to be compatible with following workflow
    exprs_df = exprs_df.transpose()

    # Transform the gene expression data (z-score by gene)
    scaled_fit = StandardScaler().fit(exprs_df)
    exprs_scaled_df_cohort = pd.DataFrame(
      scaled_fit.transform(exprs_df), index=exprs_df.index, columns=exprs_df.columns
    )
    exprs_scaled_df = exprs_scaled_df.append(exprs_scaled_df_cohort)
    print("completed scale "+cohort+" "+rna_library)

    # Shuffle input RNAseq matrix and apply classifiers
    exprs_shuffled_df_cohort = exprs_scaled_df_cohort.apply(shuffle_columns, axis=0)
    exprs_shuffled_df=exprs_shuffled_df.append(exprs_shuffled_df_cohort)
    print("completed shuffle "+cohort+" "+rna_library)


# Apply Ras Classifier

# Load RAS Classifier
file = os.path.join(
    "reference", "ras_classifier_coefficients.tsv"
)
ras_coef_df = pd.read_table(file, index_col=0)
ras_coef_df = ras_coef_df.query("abs_val > 0")

# Apply the Ras classifier to the input RNAseq matrix
ras_scores_df, ras_common_genes_df, ras_missing_genes_df = apply_classifier(
    coef_df=ras_coef_df, rnaseq_df=exprs_scaled_df
)

# Determine the extent of coefficient overlap
print(
    "There are a total of {} out of {} genes in common ({}%) between the datasets".format(
        ras_common_genes_df.shape[0],
        ras_coef_df.shape[0],
        round(ras_common_genes_df.shape[0] / ras_coef_df.shape[0] * 100, 2),
    )
)

# Which Genes are Missing?
print(ras_missing_genes_df)

# Apply Ras Classifier to Shuffled Data

# Apply the Ras classifier to the input RNAseq matrix
(
    ras_shuffle_scores_df,
    ras_shuffle_common_genes_df,
    ras_shuffle_missing_genes_df,
) = apply_classifier(coef_df=ras_coef_df, rnaseq_df=exprs_shuffled_df)

# Apply TP53 Classifier

# Load TP53 Classifier
file = os.path.join(
    "reference", "tp53_classifier_coefficients.tsv"
)
tp53_coef_df = pd.read_table(file, index_col=0)
tp53_coef_df = tp53_coef_df.query("abs_val > 0")

# Apply the TP53 classifier to the input RNAseq matrix
tp53_scores_df, tp53_common_genes_df, tp53_missing_genes_df = apply_classifier(
    coef_df=tp53_coef_df, rnaseq_df=exprs_scaled_df
)

# Determine the extent of coefficient overlap
print(
    "There are a total of {} out of {} genes in common ({}%) between the datasets".format(
        tp53_common_genes_df.shape[0],
        tp53_coef_df.shape[0],
        round(tp53_common_genes_df.shape[0] / tp53_coef_df.shape[0] * 100, 2),
    )
)

# Which Genes are Missing?
print(tp53_missing_genes_df)

# Apply TP53 Classifier to Shuffled Data

# Apply the Tp53 classifier to the input RNAseq matrix
(
    tp53_shuffle_scores_df,
    tp53_shuffle_common_genes_df,
    tp53_shuffle_missing_genes_df,
) = apply_classifier(coef_df=tp53_coef_df, rnaseq_df=exprs_shuffled_df)

# Apply NF1 Classifier

# Load NF1 Classifier
file = os.path.join(
    "reference", "nf1_classifier_coefficients.tsv"
)
nf1_coef_df = pd.read_table(file, index_col=0)
nf1_coef_df = nf1_coef_df.query("abs_val > 0")

# Apply the NF1 classifier to the input RNAseq matrix
nf1_scores_df, nf1_common_genes_df, nf1_missing_genes_df = apply_classifier(
    coef_df=nf1_coef_df, rnaseq_df=exprs_scaled_df
)

# Determine the extent of coefficient overlap
print(
    "There are a total of {} out of {} genes in common ({}%) between the datasets".format(
        nf1_common_genes_df.shape[0],
        nf1_coef_df.shape[0],
        round(nf1_common_genes_df.shape[0] / nf1_coef_df.shape[0] * 100, 2),
    )
)

# Which Genes are Missing?
print(nf1_missing_genes_df)

# Apply NF1 Classifier to Shuffled Data

# Apply the NF1 classifier to the input RNAseq matrix
(
    nf1_shuffle_scores_df,
    nf1_shuffle_common_genes_df,
    nf1_shuffle_missing_genes_df,
) = apply_classifier(coef_df=nf1_coef_df, rnaseq_df=exprs_shuffled_df)

# Combine Ras, NF1, and TP53 predictions and output to file

results_list = [
    ras_scores_df.T,
    tp53_scores_df.T,
    nf1_scores_df.T,
    ras_shuffle_scores_df.T,
    tp53_shuffle_scores_df.T,
    nf1_shuffle_scores_df.T,
]

all_results = pd.concat(results_list, axis="columns").reset_index()
all_results.columns = [
    "sample_id",
    "ras_score",
    "tp53_score",
    "nf1_score",
    "ras_shuffle",
    "tp53_shuffle",
    "nf1_shuffle",
]

filename = pd.Series([name, "classifier_scores.tsv"])
filename = filename.str.cat(sep="_")
print(filename)
file = os.path.join("results", filename)
all_results.to_csv(file, sep="\t", index=False)
