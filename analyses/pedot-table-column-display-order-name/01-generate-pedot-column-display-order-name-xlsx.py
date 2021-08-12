#!/usr/bin/env python3
import os
import pandas as pd
from utils import TSVSheet

# List of all TSV sheet names and file paths to generate the xlsx spreadsheet
tsv_sheets = [
    TSVSheet(
        "SNV gene-level",
        os.path.join("..", "snv-frequencies", "results",
                     "gene-level-snv-consensus-annotated-mut-freq.tsv")),
    TSVSheet(
        "SNV variant-level",
        os.path.join("..", "snv-frequencies", "results",
                     "variant-level-snv-consensus-annotated-mut-freq.tsv.gz")),
    TSVSheet(
        "CNV gene-level",
        os.path.join("..", "cnv-frequencies", "results",
                     "gene-level-cnv-consensus-annotated-mut-freq.tsv.gz")),
    TSVSheet(
        "Fusion gene-level",
        os.path.join("..", "fusion-frequencies", "results",
                     "putative-oncogene-fused-gene-freq.tsv.gz")),
    TSVSheet(
        "Fusion fusion-level",
        os.path.join("..", "fusion-frequencies", "results",
                     "putative-oncogene-fusion-freq.tsv.gz")),
    TSVSheet(
        "TPM stats gene-wise z-scores",
        os.path.join("..", "rna-seq-expression-summary-stats", "results",
                     "long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz")),
    TSVSheet(
        "TPM stats group-wise z-scores",
        os.path.join("..", "rna-seq-expression-summary-stats", "results",
                     "long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz"))
]
# assert all xlsx sheet names are unique
assert len(tsv_sheets) == len(set([x.xlsx_sheet_name for x in tsv_sheets]))

# Path to output xlsx spreadsheet
output_xlsx_path = os.path.join("results",
                                "pedot-table-column-display-order-name.xlsx")


def main():
    # pylint: disable=abstract-class-instantiated
    with pd.ExcelWriter(output_xlsx_path, engine="openpyxl") as xlsx_writer:
        # pylint: enable=abstract-class-instantiated
        for tsv_sheet in tsv_sheets:
            tsv_sheet.write_xlsx_sheet(xlsx_writer)


if __name__ == "__main__":
    main()
