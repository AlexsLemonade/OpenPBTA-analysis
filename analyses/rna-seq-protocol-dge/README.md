## Differential gene expression analysis comparing poly-A and ribodeplete-stranded RNA-seq protocols

**Module authors:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

### Purpose

Perform differential gene expression (DGE) analysis to compare RNA-seq libraries that are prepared using poly-A or ribodeplete-stranded protocols from the same samples.

Check whether the p-values from DGE analysis follow a uniform distribution.

Select genes that are stably expressed in poly-A and stranded RNA-seq libraries.

### Methods

1. Collapse `pbta-gene-counts-rsem-expected_count` of stranded and poly-A libraries using the OpenPBTA-analysis/analyses/collapse-rnaseq module.

### Usage

1. Change working directory to local `OpenPBTA-analysis`.
2. Download data using `bash download-data.sh`. Make sure `data/pbta-gene-counts-rsem-expected_count.polya.rds` and `data/pbta-gene-counts-rsem-expected_count.stranded.rds` are downloaded.
3. Run this analysis module in the continuous integration (CI) docker image using `./scripts/run_in_ci.sh bash analyses/rna-seq-protocol-dge/run-rna-seq-protocol-dge.sh`.

Note on downloading data:

The presented results were generated using [AlexsLemonade/OpenPBTA-analysis](https://github.com/AlexsLemonade/OpenPBTA-analysis) release-v19-20210423. The release data can be downloaded by changing the `URL`, `RELEASE` and `PREVIOUS` variables in `download-data.sh` to:

```bash
URL=${OPENPBTA_URL:-https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/data}
RELEASE=${OPENPBTA_RELEASE:-release-v19-20210423}
PREVIOUS=${OPENPBTA_RELEASE:-release-v18-20201123}
```

The [PediatricOpenTargets/OpenPBTA-analysis](https://github.com/PediatricOpenTargets/OpenPBTA-analysis) project has a separate data release, so the `download-data.sh` currently has the following variables:

```bash
URL=${OPENPBTA_URL:-https://s3.amazonaws.com/kf-openaccess-us-east-1-prd-pbta/open-targets}
RELEASE=${OPENPBTA_RELEASE:-v2}
PREVIOUS=${OPENPBTA_RELEASE:-v1}
```

The `data/pbta-gene-counts-rsem-expected_count.polya.rds` and `data/pbta-gene-counts-rsem-expected_count.stranded.rds` will be added to PediatricOpenTargets/OpenPBTA-analysis in future releases.

### Analysis scripts

#### 01-summarize_matrices.R

> this script generates the collapsed matrices as described above. In addition, this script calculates the average Pearson correlation between the values of the gene symbol that is kept and those duplicates that are discarded.

Obtained from the OpenPBTA-analysis/analyses/collapse-rnaseq module.

Input:

- `../../data/pbta-gene-counts-rsem-expected_count.polya.rds`
- `../../data/pbta-gene-counts-rsem-expected_count.stranded.rds`

Output:

- `results/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds`: collapsed poly-A read count matrix
- `results/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds`: collapsed stranded read count matrix
- `results/pbta-gene-counts-rsem-expected_count-collapsed_table.polya.rds`
- `results/pbta-gene-counts-rsem-expected_count-collapsed_table.stranded.rds`

#### 02-analyze-drops.Rmd

> this is used to display tables from 01-summarize_matrices.R

Obtained from the OpenPBTA-analysis/analyses/collapse-rnaseq module.

Input:

- `results/pbta-gene-counts-rsem-expected_count-collapsed_table.polya.rds`
- `results/pbta-gene-counts-rsem-expected_count-collapsed_table.stranded.rds`

Output:

- `results/expected_count-collapse-gene-drops.nb.html`: dropped gene report
