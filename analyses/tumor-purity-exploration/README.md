## Tumor purity exploration module

Author: SJ Spielman `@sjspielman`

This analysis modules explores the tumor purity values and potential correlates.

The full module can be run with (assuming your current directory is the root of this module):

```
bash run_tumor-purity.sh
```

Contents of this module include the following:

- `01_explore-tumor-purity.Rmd`: Notebook to explore the overall distribution of purity values and their relationship to other metadata information including diagnoses and extraction type.
This notebook outputs the file `results/rna_stranded_same-extraction.tsv`, which contains the metadata (including mapped tumor fraction) for RNA stranded samples that are from the same extraction as their DNA counterparts from which tumor fraction was measured.

- `02_explore-tumor-purity-threshold.Rmd`: Notebook to explore resulting sample sizes after filtering samples based on different thresholds for tumor purity values.
This notebook outputs the file `results/thresholded_rna_stranded_same-extraction.tsv`, which is a subsetted version of `results/rna_stranded_same-extraction.tsv` that contains only those tumors whose purity is >= the respective cancer group's median tumor purity value.