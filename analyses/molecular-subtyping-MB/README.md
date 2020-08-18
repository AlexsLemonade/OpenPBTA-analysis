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

This script first subsets the input matrices to Medulloblastoma samples only. Next, the two matrices containing MB samples are merged together using common genes. Then, if the value of the input parameter `batch_col` is set, a batch correction is performed using that column from the clinical file. Finally, the script runs either the [MM2S](https://github.com/cran/MM2S) package or the medulloblastoma classifier as described [here](https://github.com/d3b-center/medullo-classifier-package) depending on the `--method` parameter.

3. Output: 

```
results/mb-molecular-subtypes*.rds
```

The results in the rds object contain the samples, best.fit which is medulloblastoma subtype assigned to the sample and associated p-value.

#### 02-compare-classes.R 

1. Input

```
results/mb-molecular-subtypes*.rds
```

2. Function:

This script takes in the expected classification from pathology reports and merges it with the observed classification obtained after running the classifier.

3. Output

```
results/comparison-*.rds
```

### Running the analysis

```sh 
bash run-molecular-subtyping-MB*.sh
```



