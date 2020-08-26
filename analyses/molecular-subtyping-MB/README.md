## Molecular Subtype Classification (MB)

**Module authors:** Komal S. Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to leverage the R packages [medulloPackage](https://github.com/d3b-center/medullo-classifier-package) and [MM2S](https://github.com/cran/MM2S) that utilize expression data from RNA-seq or array to classify the medulloblastoma (MB) samples into four subtypes i.e Group3, Group4, SHH, WNT. The input for both classifiers is a log-normalized FPKM matrix with gene symbols as rownames in case of medulloPackage and entrez ids as rownames in case of MM2S.  

We use the polyA selected (n = 1) and rRNA depleted (n = 121) MB samples as input. Because the classifiers do not work on single-sample, we merge the two matrices by common gene symbols and use the merged matrix as input. To see if batch correction of library type as any effect on the classification, we also use a second input matrix that has been batch corrected for library type using the sva package.

### Analysis scripts

#### 00-filter-and-batch-correction.R

1. Inputs from data download

```
pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
pbta-histologies.tsv
```

2. Function

This script first subsets the input expression matrices to MB samples only. Next, the poly-A and rRNA depleted input matrices containing only MB samples are merged together using common genes. Using this merged matrix, we then generate two input matrices for molecular subtype classification: 

	1. uncorrected log-normalized FPKM matrix and 
	2. batch corrected, log-normalized FPKM matrix. 

The batch correction uses the column from the clinical file that matches the value of `--batch_col` as the batch variable and computes a batch correction using the R package sva.

3. Output: 

```
# uncorrected matrix
results/medulloblastoma-exprs.rds

# batch corrected matrix
results/medulloblastoma-exprs-batch-corrected.rds
```

#### 01-classify-mb.R 

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

#### 02-compare-classes.Rmd

1. Input

```
# expected output from pathology reports
input/expected_class.rds

# observed output from 01-classify-mb.R 
results/mb-classified.rds
``` 

2. Function:

This notebook summarizes the performance of the two classifiers on batch corrected and uncorrected expression matrix obtained after running 01-classify-mb.R. 

The observed subtypes obtained from each classifier are compared to the expected subtypes in the pathology report in order to determine the classifier accuracy. % Accuracy is calculated by matching observed and expected subtypes only where expected subtype information is available. In case of ambiguous subtypes, a match is determined only if the observed subtype matches with any one of the expected subtypes. 

The pathology report has subtype information on 32/122 (26.2%) samples. Following is the breakdown of pathology identified subtypes:

| pathology_subtype | freq |
|-------------------|------|
| Group 3 or 4      | 11   |
| Group 4           | 2    |
| non-WNT           | 1    |
| SHH               | 12   |
| WNT               | 6    |

For each input expression matrix (i.e. uncorrected and batch-corrected), the molecular subtype is determined by taking a consensus of the two classifiers. If the classifiers do not agree, the sample is treated as unclassified.

3. Output

```
# html output
02-compare-classes.html
```

### Running the full analysis

```sh 
bash run-molecular-subtyping-mb.sh
```



