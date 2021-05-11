## Molecular Subtype Classification (MB)

**Module authors:** Komal S. Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to leverage the R packages [medulloPackage](https://github.com/d3b-center/medullo-classifier-package) and [MM2S](https://github.com/cran/MM2S) that utilize expression data from RNA-seq or array to classify the medulloblastoma (MB) samples into four subtypes i.e Group3, Group4, SHH, WNT. The input for both classifiers is a log-normalized FPKM matrix with gene symbols as rownames in case of medulloPackage and entrez ids as rownames in case of MM2S.  

We use the polyA selected (n = 1) and rRNA depleted (n = 121) MB samples as input. Because the classifiers do not work on single-sample, we merge the two matrices by common gene symbols and use the merged matrix as input. To see if batch correction of library type as any effect on the classification, we also use a second input matrix that has been batch corrected for library type using the sva package.

### Running the full analysis

This runs 01 - 03 scripts to create all the output in the `results/` folder.

```sh
bash run-molecular-subtyping-mb.sh
```

### Analysis scripts

#### 00-mb-select-pathology-dx.Rmd

1. Inputs

```
pbta-histologies.tsv
```

2. Function

_This 00 script is not run by calling run-molecular-subtyping-mb.sh_ but can be run separately using this command:

```sh
Rscript -e "rmarkdown::render('analyses/molecular-subtyping-MB/00-mb-select-pathology-dx.Rmd', clean = TRUE)"
```

This Rmd checks the alignment of `Medulloblastoma` labels across fields:  `pathology_diagnosis`, `integrated_diagnosis`, and `short_histology`.
It is tied to a specific release (`release-v18-20201123`).
This creates a terms JSON which is used in the other scripts for subsetting.

3. Output:

A medulloblastoma terms JSON file:
`molecular-subtyping-MB/inputs/mb_subtyping_path_dx_strings.json`

#### 01-filter-and-batch-correction.R

1. Inputs

```
# the rna-seq expression files
data/pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds

# medulloblastoma terms file
molecular-subtyping-MB/inputs/mb_subtyping_path_dx_strings.json
```

2. Function

This script first subsets the input expression matrices to MB samples only. Next, the poly-A and rRNA depleted input matrices containing only MB samples are merged together using common genes. Using this merged matrix, we then generate two input matrices for molecular subtype classification:

	1. uncorrected log-normalized FPKM matrix and
	2. batch corrected, log-normalized FPKM matrix.

The batch correction uses the column from the clinical file that matches the value of `--batch_col` as the batch variable and computes a batch correction using the R package sva.

3. Output:

```
# subset clinical file to medulloblastoma biospecimens only
input/subset-mb-clinical.tsv

# uncorrected matrix
results/medulloblastoma-exprs.rds

# batch corrected matrix
results/medulloblastoma-exprs-batch-corrected.rds
```

#### 02-classify-mb.R

1. Input

```
# uncorrected matrix
results/medulloblastoma-exprs.rds

# batch corrected matrix
results/medulloblastoma-exprs-batch-corrected.rds
```

2. Function:

This script runs the two classifiers on both uncorrected and batch-corrected input matrices. In order to run MM2S, the script utilizes the R package `org.Hs.eg.db`  to convert gene symbols to Entrez ids.

3. Output

```
results/mb-classified.rds
```

The .rds object contains a list of dataframes with outputs corresponding to the four runs as described above. Each dataframe contains 5 columns: sample (Kids_First_Biospecimen_ID), best.fit (i.e. medulloblastoma subtype assigned to the sample), classifier (MM2S or medullo-classifier), dataset (corrected or uncorrected matrix) and score (in case of MM2S) or p-value (in case of medulloPackage).  

#### 03-compare-classes.Rmd

1. Input

```
# subset clinical file to medulloblastoma biospecimens only
input/subset-mb-clinical.tsv

# expected output from pathology reports
input/openPBTA-mb-pathology-subtypes.rds

# observed output from 01-classify-mb.R
results/mb-classified.rds
```

TODO: `input/openPBTA-mb-pathology-subtypes.rds` doesn't have much information surrounding it, but [this is issue has been tracked](https://github.com/AlexsLemonade/OpenPBTA-analysis/issues/746) and this README can be updated once that has been addressed.

2. Function:

This notebook summarizes the performance of the two classifiers on batch corrected and uncorrected expression matrix obtained after running 01-classify-mb.R.

The observed subtypes obtained from each classifier are compared to the expected subtypes in the pathology report in order to determine the classifier accuracy. % Accuracy is calculated by matching observed and expected subtypes only where expected subtype information is available. In case of ambiguous subtypes, a match is determined only if the observed subtype matches with any one of the expected subtypes.

The pathology report has subtype information on 43/122 (35.2%) samples. Following is the breakdown of pathology identified subtypes:

| pathology_subtype | freq |
|-------------------|------|
| Group 3 or 4      | 14   |
| Group 4           | 5    |
| non-WNT           | 9    |
| SHH               | 10   |
| WNT               | 5    |

For each input expression matrix (i.e. uncorrected and batch-corrected), the molecular subtype is determined by taking a consensus of the two classifiers. For cases in which the classifiers do not agree, non-ambiguous subtypes from `pathology_subtype` field is taken and the remaining samples are treated as unclassified.

For each sample_id with multiple RNA samples, the following logic is implemented:
- IF there are two RNA specimens from the same event (`sample_id`) with the same `tumor_descriptor`, and
- IF one of the two RNA specimens has a consensus match from the MB classifiers, and
- IF the RNA specimen with mismatched classifications includes a match to the subtype of the consensus of the other RNA specimen,
- THEN propagate the subtype to the second RNA specimen.

For the consensus corrected output, 37/43 (86.05%) and for the consensus uncorrected output, 38/43 (88.37%) match with the reported subtype. The following table shows the difference between the two outputs:

BS_V96WVE3Z is correctly predicted by consensus uncorrected output but not by batch-corrected output:

| Consensus Output | Kids_First_Biospecimen_ID | pathology_subtype | MM2S_best_fit | medulloPackage_best_fit | molecular_subtype | match |
|------------------|---------------------------|-------------------|---------------|-------------------------|-------------------|-------|
| Uncorrected      | BS_V96WVE3Z               | SHH               | SHH           | SHH                     | SHH               | TRUE  |
| Corrected        | BS_V96WVE3Z               | SHH               | Group3        | SHH                     | NA                | NA    |

3. Output

The markdown produces one html notebook and a tab-delimited file containing RNA and DNA identifiers mapped to the consensus molecular subtype obtained for each input matrix.

```
# html output
02-compare-classes.html

# tsv files containing RNA and DNA identifiers mapped to molecular_subtype
# uncorrected input consensus output
results/MB_molecular_subtype.tsv

# batch corrected input consensus output
results/MB_batchcorrected_molecular_subtype.tsv
```
