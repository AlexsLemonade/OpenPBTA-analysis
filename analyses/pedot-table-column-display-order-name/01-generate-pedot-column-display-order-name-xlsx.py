#!/usr/bin/env python3
import os
import datetime
import pandas as pd
from utils import TSVSheet
from utils import update_xlsx_datetime

# List of all TSV file paths to generate the xlsx spreadsheet
tsv_paths = [
    os.path.join("..", "snv-frequencies", "results",
                 "gene-level-snv-consensus-annotated-mut-freq.tsv.gz"),
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
    # date time tuple for xlsx file creation
    xlsx_creation_datetime_tuple = (2021, 8, 12, 17, 20, 7)
    with pd.ExcelWriter(output_xlsx_path, engine="openpyxl") as xlsx_writer:
        # Set creation time to make output file identically reproducible
        xlsx_creation_datetime = datetime.datetime(
            *xlsx_creation_datetime_tuple)
        xlsx_writer.book.properties.created = xlsx_creation_datetime
        xlsx_writer.book.properties.modified = xlsx_creation_datetime
        for tp in tsv_paths:
            tsv_sheet = TSVSheet(tp)
            tsv_sheet.write_xlsx_sheet(xlsx_writer)
    # xlsx file is actually a zip file
    # Set datetimes in xlsx zip file to make xlsx file identically
    # reproducible
    update_xlsx_datetime(output_xlsx_path, xlsx_creation_datetime_tuple)


if __name__ == "__main__":
    main()
