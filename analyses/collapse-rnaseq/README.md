## Collapse RNA-seq matrices

Many downstream modules require RNA-seq data that have gene symbols as gene identifiers.
Multiple Ensembl gene identifiers map to the same gene symbol in a small proportion of cases.
This module contains the script used to create RSEM matrices with duplicate gene symbols removed, by first filtering to only genes with [FPKM/TPM] > 0 and then selecting the instance of the gene symbol with the maximum mean [FPKM/TPM/Expected_count] value across samples.
It produces the files with below naming patterns included in the data download:

```
gene-[expression/counts]-rsem-[fpkm/tpm/expected_count]-collapsed.rds
```

This script also runs the notebook for analysis of the dropped genes, and produces a markdown report.

```
gene-[expression/counts]-rsem-[fpkm/tpm/expected_count]-collapsed_table.rds
```

To run the steps that generate the matrices and display the results of a correlation analysis, use the following command (assuming you are in this directory):

```sh
bash run-collapse-rnaseq.sh -q <quatification_type> -x <'type_of_data'>
```

Acceptable values for the options are below:

<quatification_type> : 'fpkm', 'tpm', 'expected_count' <br>
<type_of_data> : 'expression', 'counts' <br>


### R scripts and notebook

* `00-create-rsem-files.R` - this script is used to split RSEM files into two matrices, one for each library strategy. This step occurs upstream of this repository to produce the RSEM files included in the data download. 
This script is not run via `run-collapse-rnaseq.sh`.
* `01-summarize_matrices.R` - this script generates the collapsed matrices as described above.
In addition, this script calculates the average Pearson correlation between the values of the gene symbol that is kept and those duplicates that are discarded.
* `02-analyze-drops.Rmd` - this is used to display tables from `01-summarize_matrices.R`.
