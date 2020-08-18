## Molecular Subtype Classification (MB)

**Module authors:** Komal S. Rathi ([@komalsrathi](https://github.com/komalsrathi))

### Description

The goal of this analysis is to use the R package [medulloPackage](https://github.com/d3b-center/medullo-classifier-package) to leverage expression data from RNA-seq or array to classify the medulloblastoma samples into four subtypes i.e Group3, Group4, SHH, WNT.

### Analysis scripts

#### 01-classify-mb.R

1. Inputs from data download

```
pbta-gene-expression-rsem-fpkm-collapsed.polya.rds
pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds
pbta-histologies.tsv
```

2. Function

This script first subsets the input matrices to Medulloblastoma samples only. Next, the two matrices containing MB samples are merged together using common genes. Then, if the value of the input parameter `batch_col` is set, a batch correction is performed using that column from the clinical file. Finally, the script runs the medulloblastoma classifier as described [here](https://github.com/d3b-center/medullo-classifier-package).

3. Output: 

```
# with batch correction
results/mb-molecular-subtypes-v1.rds

# no batch correction
results/mb-molecular-subtypes-v2.rds
```

The results in the rds object contain the samples, best.fit which is medulloblastoma subtype assigned to the sample and associated p-value.

#### 02-compare-classes.R 

1. Input

```
# using batch correction
results/mb-molecular-subtypes-v1.rds

# no batch correction
results/mb-molecular-subtypes-v2.rds
```

2. Function:

This script takes in the expected classification from pathology reports and merges it with the observed classification obtained after running the classifier.

3. Output

```
# using batch correction
results/comparison_v1.rds

# no batch correction
results/comparison_v2.rds
```

### Running the analysis

```sh
# with batch correction using library type as batch 
bash run-molecular-subtyping-MB-v1.sh

# without batch correction
bash run-molecular-subtyping-MB-v2.sh
```



