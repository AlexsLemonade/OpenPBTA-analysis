## Molecular Subtype Classification (MB)

**Module authors:** Komal S. Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to leverage the R package [medulloPackage](https://github.com/d3b-center/medullo-classifier-package) that utilizes expression data from RNA-seq or array to classify the medulloblastoma (MB) samples into four subtypes i.e Group3, Group4, SHH, WNT. The input for medulloPackage classifier is a log-normalized FPKM matrix with gene symbols as rownames.  
We use the polyA selected (n = 1) and rRNA depleted (n = 121) MB samples as input. To see if batch correction of library type as any effect on the classification, we also use a second input matrix that has been batch corrected for library type using the sva package.

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
This creates a terms JSON which is used in the other scripts for subsetting.

3. Output:

A medulloblastoma terms JSON file:
`molecular-subtyping-MB/inputs/mb_subtyping_path_dx_strings.json`

#### 01-filter-and-batch-correction.R

1. Inputs

```
# the rna-seq expression files
data/gene-expression-rsem-tpm-collapsed.rds

# histologies file
data/histologies.tsv

# medulloblastoma terms file
molecular-subtyping-MB/inputs/mb_subtyping_path_dx_strings.json
```

2. Function

This script first subsets the input gene expression matrix to MB samples only. Next, it generates two input matrices for molecular subtype classification:

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

This script runs the medulloPackage classifier on both uncorrected and batch-corrected input matrices.

3. Output

```
results/mb-classified.rds
```

The .rds object contains a list of dataframes with outputs corresponding to the two runs as described above. Each dataframe contains 5 columns: `sample` (Kids_First_Biospecimen_ID), `best.fit` (i.e. medulloblastoma subtype assigned to the sample), `classifier` (medullo-classifier), `dataset` (corrected or uncorrected matrix) and `p.value` (in case of medulloPackage).  

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

For each input expression matrix (i.e. uncorrected and batch-corrected), the molecular subtype is determined by using the prediction obtained from `medulloPackage` classifier. 

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
