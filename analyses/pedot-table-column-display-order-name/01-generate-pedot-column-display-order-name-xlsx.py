#!/usr/bin/env python3
import os
from utils import TSV
import pandas as pd

# List of all TSV file paths to generate the xlsx spreadsheet
tsv_paths = [
    os.path.join("..", "snv-frequencies", "results",
                 "gene-level-snv-consensus-annotated-mut-freq.tsv"),
    os.path.join("..", "snv-frequencies", "results",
                 "variant-level-snv-consensus-annotated-mut-freq.tsv.gz"),
    os.path.join("..", "cnv-frequencies", "results",
                 "gene-level-cnv-consensus-annotated-mut-freq.tsv.gz"),
    os.path.join("..", "fusion-frequencies", "results",
                 "putative-oncogene-fused-gene-freq.tsv.gz"),
    os.path.join("..", "fusion-frequencies", "results",
                 "putative-oncogene-fusion-freq.tsv.gz"),
    os.path.join("..", "rna-seq-expression-summary-stats", "results",
                 "long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz"),
    os.path.join("..", "rna-seq-expression-summary-stats", "results",
                 "long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz")
]

# Path to output xlsx spreadsheet
output_xlsx_path = os.path.join("results",
                                "pedot-table-column-display-order-name.xlsx")


def main():
    with pd.ExcelWriter(output_xlsx_path, engine="openpyxl") as xlsx_writer:
        for path in tsv_paths:
            tsv = TSV(path)
            tsv.write_xlsx_sheet(xlsx_writer)


if __name__ == "__main__":
    main()
