## Tumor purity exploration module

Author: SJ Spielman `@sjspielman`

This analysis modules explores the tumor purity values and potential correlates.

The full module can be run with `bash run_tumor-purity`.

Contents include the following:

- `01-explore-tumor-purity.Rmd`: Notebook to explore the overall distribution of purity values and their relationship to other metadata information including diagnoses and extraction type.
This notebook outputs two files:
    -`results/metadata_rna_stranded_same-extraction.tsv`, which contains the metadata (including mapped tumor fraction) for RNA stranded samples that are from the same extraction as their DNA counterparts from which tumor fraction was measured.
    -`results/metadata_thresholded_rna_stranded_same-extraction.tsv`, is a subsetted version of the above TSV that contains only those tumors whose purity is >= the respective cancer group's median tumor purity value.