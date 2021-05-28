## Differential gene expression analysis comparing poly-A and ribodeplete-stranded RNA-seq protocols

**Module authors:** Yuanchao Zhang ([@logstar](https://github.com/logstar))

### Purpose

Perform differential gene expression (DGE) analysis to compare RNA-seq libraries that are prepared using poly-A or ribodeplete-stranded protocols from the same samples.

Check whether the p-values from DGE analysis follow a uniform distribution.

Select genes that are stably expressed in poly-A and stranded RNA-seq libraries.

### Methods

1. Collapse `pbta-gene-counts-rsem-expected_count` of stranded and poly-A libraries using the OpenPBTA-analysis/analyses/collapse-rnaseq module.
2. Select `Kids_First_Biospecimen`s with the same `sample_id`s and both stranded and poly-A `RNA_library`s.
    Index | Kids_First_Biospecimen_ID | sample_id | experimental_strategy | RNA_library | cohort
    ------|---------------------------|-----------|-----------------------|-------------|--------
    1     | BS_HE0WJRW6               | 7316-1455 | RNA-Seq               | stranded    | CBTN
    2     | BS_HWGWYCY7               | 7316-1455 | RNA-Seq               | poly-A      | CBTN
    3     | BS_SHJA4MR0               | 7316-161  | RNA-Seq               | stranded    | CBTN
    4     | BS_X0XXN9BK               | 7316-161  | RNA-Seq               | poly-A      | CBTN
    5     | BS_FN07P04C               | 7316-255  | RNA-Seq               | stranded    | CBTN
    6     | BS_W4H1D4Y6               | 7316-255  | RNA-Seq               | poly-A      | CBTN
    7     | BS_8QB4S4VA               | 7316-536  | RNA-Seq               | stranded    | CBTN
    8     | BS_QKT3TJVK               | 7316-536  | RNA-Seq               | poly-A      | CBTN
    9     | BS_7WM3MNZ0               | A16915    | RNA-Seq               | poly-A      | PNOC003
    10    | BS_KABQQA0T               | A16915    | RNA-Seq               | stranded    | PNOC003
    11    | BS_68KX6A42               | A18777    | RNA-Seq               | poly-A      | PNOC003
    12    | BS_D7XRFE0R               | A18777    | RNA-Seq               | stranded    | PNOC003
3. Run edgeR `exactTest`, edgeR `glmLRT`, edgeR `glmQLFTest`, and DESeq2 `nbinomWaldTest` comparing stranded and poly-A RNA-seq `rsem-expected_count`s, using an adapted workflow from the [OMPARE project](https://github.com/d3b-center/OMPARE/blob/master/code/patient_level_analyses/utils/rnaseq_edger_normalizations.R#L37-L41). RSEM expected read counts are normalized using either UQ-pgQ2 or TMM. The UQ-pgQ2 normalization method has lower false positive rate in DGE testing (<https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6502-7>).
4. Select genes that are stably expressed in poly-A and stranded RNA-seq libraries by the following criteria, and the filtered output `results/stranded_vs_polya_stably_exp_genes.csv` is ordered by the coefficient of variations.
   - `coefficient of variation < quantile(coefficient of variation, 0.25)`
   - `exactTest`, `glmLRT`, `glmQLFTest`, and `nbinomWaldTest` p-values all > 0.05
   - `abs(logFC) < median(abs(logFC)`
   - `quantile(average_logCPM, 0.25) < average_logCPM < quantile(average_logCPM, 0.75)`
5. Run DESeq2 standard DGE analysis on RSEM expected read counts comparing stranded and poly-A RNA-seq `rsem-expected_count`s.

### Results

#### DGE QL F-test

![up_qlft_pvals](plots/uqpgq2_normalized/stranded_vs_polya_dge_ql_ftest_pvals_histogram.png)

![tmm_qlft_pvals](plots/tmm_normalized/stranded_vs_polya_dge_ql_ftest_pvals_histogram.png)

DGE QL F-test result tables:

- UQ-pgQ2 normalized: `results/uqpgq2_normalized/stranded_vs_polya_dge_ql_ftest_res.csv`
- TMM normalized: `results/tmm_normalized/stranded_vs_polya_dge_ql_ftest_res.csv`

#### DGE LRT

![up_lrt_pvals](plots/uqpgq2_normalized/stranded_vs_polya_dge_lrt_pvals_histogram.png)

![tmm_lrt_pvals](plots/tmm_normalized/stranded_vs_polya_dge_lrt_pvals_histogram.png)

DGE LRT result tables:

- UQ-pgQ2 normalized: `results/uqpgq2_normalized/stranded_vs_polya_dge_lrt_res.csv`
- TMM normalized: `results/tmm_normalized/stranded_vs_polya_dge_lrt_res.csv`

#### DGE exactTest

![up_exact_test_pvals](plots/uqpgq2_normalized/stranded_vs_polya_dge_exact_test_pvals_histogram.png)

![tmm_exact_test_pvals](plots/tmm_normalized/stranded_vs_polya_dge_exact_test_pvals_histogram.png)

DGE exactTest result tables:

- UQ-pgQ2 normalized: `results/uqpgq2_normalized/stranded_vs_polya_dge_exact_test_res.csv`
- TMM normalized: `results/tmm_normalized/stranded_vs_polya_dge_exact_test_res.csv`

#### DGE DESeq2 nbinomWaldTest

![up_exact_test_pvals](plots/uqpgq2_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png)

![tmm_exact_test_pvals](plots/tmm_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png)

![ds_exact_test_pvals](plots/deseq2_rle_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png)

DGE exactTest result tables:

- UQ-pgQ2 normalized: `results/uqpgq2_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv`
- TMM normalized: `results/tmm_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv`
- DESeq2 standard DGE analysis: `results/deseq2_rle_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv`

#### Stably expressed genes

##### UQ-pgQ2 normalization

The table of stably expressed genes are saved at `results/uqpgq2_normalized/stranded_vs_polya_stably_exp_genes.csv`. The biospecimen ID columns are normalized counts.

Top 30 stably expressed gene boxplots are at `plots/uqpgq2_normalized/stably_exp_gene_protocol_diff_boxplot`. Top 3 stably expressed gene boxplots:

![up_seg_bp1](plots/uqpgq2_normalized/stably_exp_gene_protocol_diff_boxplot/RAB5B_normalized_count_protocol_diff_boxplot.png)

![up_seg_bp2](plots/uqpgq2_normalized/stably_exp_gene_protocol_diff_boxplot/GORASP1_normalized_count_protocol_diff_boxplot.png)

![up_seg_bp3](plots/uqpgq2_normalized/stably_exp_gene_protocol_diff_boxplot/SENP2_normalized_count_protocol_diff_boxplot.png)

Histograms used for choosing stably expressed genes:

![up_norm_cnt_cv_hist](plots/uqpgq2_normalized/normalized_count_gene_cv_histogram.png)

![up_logfc_hist](plots/uqpgq2_normalized/edger_logfc_histogram.png)

##### TMM normalization

The table of stably expressed genes are saved at `results/tmm_normalized/stranded_vs_polya_stably_exp_genes.csv`. The biospecimen ID columns are normalized counts.

Top 30 stably expressed gene boxplots are at `plots/tmm_normalized/stably_exp_gene_protocol_diff_boxplot`. Top 3 stably expressed gene boxplots:

![tmm_seg_bp1](plots/tmm_normalized/stably_exp_gene_protocol_diff_boxplot/GORASP1_normalized_count_protocol_diff_boxplot.png)

![tmm_seg_bp2](plots/tmm_normalized/stably_exp_gene_protocol_diff_boxplot/RPUSD1_normalized_count_protocol_diff_boxplot.png)

![tmm_seg_bp3](plots/tmm_normalized/stably_exp_gene_protocol_diff_boxplot/NUP58_normalized_count_protocol_diff_boxplot.png)

Histograms used for choosing stably expressed genes:

![tmm_norm_cnt_cv_hist](plots/tmm_normalized/normalized_count_gene_cv_histogram.png)

![tmm_logfc_hist](plots/tmm_normalized/edger_logfc_histogram.png)

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

### Module structure

```text
.
├── 01-summarize_matrices.R
├── 02-analyze-drops.Rmd
├── 03-protocol-dge-seg.R
├── 04-deseq2-protocol-dge.R
├── README.md
├── plots
│   ├── deseq2_rle_normalized
│   │   └── stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png
│   ├── tmm_normalized
│   │   ├── edger_logfc_histogram.png
│   │   ├── estimated_dispersions.png
│   │   ├── glm_ql_fit_deviances.png
│   │   ├── normalized_count_gene_cv_histogram.png
│   │   ├── stably_exp_gene_protocol_diff_boxplot/*.png
│   │   ├── stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png
│   │   ├── stranded_vs_polya_dge_exact_test_pvals_histogram.png
│   │   ├── stranded_vs_polya_dge_lrt_pvals_histogram.png
│   │   └── stranded_vs_polya_dge_ql_ftest_pvals_histogram.png
│   └── uqpgq2_normalized
│       ├── edger_logfc_histogram.png
│       ├── estimated_dispersions.png
│       ├── glm_ql_fit_deviances.png
│       ├── normalized_count_gene_cv_histogram.png
│       ├── stably_exp_gene_protocol_diff_boxplot/*.png
│       ├── stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png
│       ├── stranded_vs_polya_dge_exact_test_pvals_histogram.png
│       ├── stranded_vs_polya_dge_lrt_pvals_histogram.png
│       └── stranded_vs_polya_dge_ql_ftest_pvals_histogram.png
├── results
│   ├── deseq2_rle_normalized
│   │   └── stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv
│   ├── expected_count-collapse-gene-drops.nb.html
│   ├── pbta-gene-counts-rsem-expected_count-collapsed.polya.rds
│   ├── pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds
│   ├── pbta-gene-counts-rsem-expected_count-collapsed_table.polya.rds
│   ├── pbta-gene-counts-rsem-expected_count-collapsed_table.stranded.rds
│   ├── tmm_normalized
│   │   ├── stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv
│   │   ├── stranded_vs_polya_dge_exact_test_res.csv
│   │   ├── stranded_vs_polya_dge_lrt_res.csv
│   │   ├── stranded_vs_polya_dge_ql_ftest_res.csv
│   │   └── stranded_vs_polya_stably_exp_genes.csv
│   └── uqpgq2_normalized
│       ├── stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv
│       ├── stranded_vs_polya_dge_exact_test_res.csv
│       ├── stranded_vs_polya_dge_lrt_res.csv
│       ├── stranded_vs_polya_dge_ql_ftest_res.csv
│       └── stranded_vs_polya_stably_exp_genes.csv
└── run-rna-seq-protocol-dge.sh
```

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

#### 03-protocol-dge-seg.R

This script performs edgeR DGE analysis to compare RNA-seq libraries that are prepared using poly-A or ribodeplete-stranded protocols from the same samples. This script also selects genes that are stably expressed in poly-A and ribodeplete-stranded RNA-seq libraries. The analysis code is adapted from the [OMPARE project](https://github.com/d3b-center/OMPARE/blob/master/code/patient_level_analyses/utils/rnaseq_edger_normalizations.R#L37-L41).

Example usage:

```bash
Rscript --vanilla '03-protocol-dge-seg.R' -n 'uqpgq2'
```

Input:

- `results/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds`: collapsed RSEM expected count matrix of poly-A RNA-seq libraries generated by `01-summarize_matrices.R`.
- `results/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds`: collapsed RSEM expected count matrix of stranded RNA-seq libraries generated by `01-summarize_matrices.R`.

Parameters:

- normalization method: tmm or uqpgq2.

Output:

- `results/NORMALIZATION_METHOD/stranded_vs_polya_dge_exact_test_res.csv`: edgeR exactTest result table
- `results/NORMALIZATION_METHOD/stranded_vs_polya_dge_lrt_res.csv`: edgeR LRT result table
- `results/NORMALIZATION_METHOD/stranded_vs_polya_dge_ql_ftest_res.csv`: edgeR QL F-test result table
- `results/NORMALIZATION_METHOD/stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv`: DESeq2 nbinomWaldTest result table
- `results/NORMALIZATION_METHOD/stranded_vs_polya_stably_exp_genes.csv`: stably expressed gene table
- `plots/NORMALIZATION_METHOD/stranded_vs_polya_dge_exact_test_pvals_histogram.png`: histogram of edgeR exactTest p-values
- `plots/NORMALIZATION_METHOD/stranded_vs_polya_dge_lrt_pvals_histogram.png`: histogram of edgeR LRT p-values
- `plots/NORMALIZATION_METHOD/stranded_vs_polya_dge_ql_ftest_pvals_histogram.png`: histogram of edgeR QL F-test p-values
- `plots/NORMALIZATION_METHOD/stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png`: histogram of DESeq2 nbinomWaldTest p-values
- `plots/NORMALIZATION_METHOD/edger_logfc_histogram.png`: histogram of edgeR computed log fold changes
- `plots/NORMALIZATION_METHOD/normalized_count_gene_cv_histogram.png`: histogram of coefficient of variations of normalized read counts
- `plots/NORMALIZATION_METHOD/stably_exp_gene_protocol_diff_boxplot/*.png`: normalized read count boxplots

#### 04-edger-protocol-dge.R

This script performs DESeq2 standard DGE analysis to compare RNA-seq libraries that are prepared using poly-A or ribodeplete-stranded protocols from the same samples.

Usage:

```bash
Rscript --vanilla '04-deseq2-protocol-dge.R'
```

Input:

- `results/pbta-gene-counts-rsem-expected_count-collapsed.polya.rds`: collapsed RSEM expected count matrix of poly-A RNA-seq libraries generated by `01-summarize_matrices.R`.
- `results/pbta-gene-counts-rsem-expected_count-collapsed.stranded.rds`: collapsed RSEM expected count matrix of stranded RNA-seq libraries generated by `01-summarize_matrices.R`.

Output:

- `results/deseq2_rle_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_res.csv`: DESeq2 nbinomWaldTest result table
- `plots/deseq2_rle_normalized/stranded_vs_polya_dge_deseq2_nbinom_wald_test_pvals_histogram.png`: histogram of DESeq2 nbinomWaldTest p-values
