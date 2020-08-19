## Molecular Subtype Classification (MB)

**Module authors:** Komal S. Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to utilize the R packages [medulloPackage](https://github.com/d3b-center/medullo-classifier-package) and [MM2S](https://github.com/cran/MM2S) that leverage expression data from RNA-seq or array to classify the medulloblastoma samples into four subtypes i.e Group3, Group4, SHH, WNT.

### Analysis scripts

#### 01-classify-mb.R

1. Inputs from data download

```
pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
pbta-histologies.tsv
```

2. Function

This script first subsets the input expression matrices to Medulloblastoma (MB) samples only. Next, the poly-A and rRNA depleted input matrices containing only MB samples are merged together using common genes. If the value of the input parameter `--batch_col` is set, the script uses the column that matches the value of `--batch_col` from the clinical file as a batch and computes a batch correction using the R package sva. Finally, the script runs either the [MM2S](https://github.com/cran/MM2S) package or the medulloblastoma classifier as described [here](https://github.com/d3b-center/medullo-classifier-package) depending on the `--method` parameter.

3. Output: 

```
results/mb-molecular-subtypes*.rds
```

The .rds object contains a dataframe with samples, best.fit (i.e. medulloblastoma subtype assigned to the sample) and associated p-value (in case of medullo-classifier-package) or score (in case of MM2S).

#### 02-compare-classes.R 

1. Input

```
results/mb-molecular-subtypes*.rds
```

2. Function:

This script takes in the expected classification from pathology reports and merges it with the observed classification obtained after running 01-classify-mb.R.

3. Output

```
results/comparison-*.rds
```

### Running the analysis

```sh 
# classify using MM2S
bash run-molecular-subtyping-MB-MM2S.sh
# classify using MM2S with batch correction on RNA_library
bash run-molecular-subtyping-MB-batch-correct-MM2S.sh

# classify using medullo-classifier
bash run-molecular-subtyping-MB-medullo-classifier.sh
# classify using medullo-classifier with batch correction on RNA_library
bash run-molecular-subtyping-MB-batch-correct-medullo-classifier.sh
```



