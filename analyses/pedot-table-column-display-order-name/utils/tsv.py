import os
from typing import List

import pandas as pd

# Variant-level SNV TSV file name, which is used for implementing the following
# known column display order change:
# - Move Variant_ID_hg38 and Protein_change columns after the
#   Gene_symbol column
VARIANT_SNV_FN = "variant-level-snv-consensus-annotated-mut-freq.tsv.gz"
# TSV filename to xlsx sheet name look up table
TSV_FN_XLSX_SN_LUT = {
    "gene-level-snv-consensus-annotated-mut-freq.tsv":
        "SNV gene-level",
    VARIANT_SNV_FN:
        "SNV variant-level",
    "gene-level-cnv-consensus-annotated-mut-freq.tsv.gz":
        "CNV gene-level",
    "putative-oncogene-fused-gene-freq.tsv.gz":
        "Fusion gene-level",
    "putative-oncogene-fusion-freq.tsv.gz":
        "Fusion fusion-level",
    "long_n_tpm_mean_sd_quantile_gene_wise_zscore.tsv.gz":
        "TPM stats gene-wise z-scores",
    "long_n_tpm_mean_sd_quantile_group_wise_zscore.tsv.gz":
        "TPM stats group-wise z-scores"
}

# assert all xlsx sheet names are unique
assert len(TSV_FN_XLSX_SN_LUT) == len(set(TSV_FN_XLSX_SN_LUT.values()))


class TSVSheet:
    """TSV sheet for generating PedOT column display order and name xlsx sheet.

    Attributes:
        xlsx_sheet_name: the name of the xlsx sheet as a string.
        col_disp_order_name_df: the column display order and name data frame as
            a pd.DataFrame.

    Raises:
        ValueError: tsv_filepath is invalid or has no corresponding .jsonl.gz
            file.
    """
    def __init__(self, tsv_filepath: str) -> None:
        """Initializes TSV with a TSV file path."""
        self._tsv_filepath = tsv_filepath
        self._tsv_basename = os.path.basename(tsv_filepath)
        self.xlsx_sheet_name = self._get_xlsx_sheet_name()
        self._jsonl_filepath = self._get_jsonl_filepath()
        self._df = self._read_tsv()
        self._col_names = self._df.columns.tolist()
        self.col_disp_order_name_df = self._get_col_disp_order_name_df()

    def _get_xlsx_sheet_name(self) -> str:
        """Get xlsx_sheet_name attribute"""
        if self._tsv_basename not in TSV_FN_XLSX_SN_LUT:
            raise ValueError("Unknown TSV file name {}".format(
                self._tsv_basename))
        return TSV_FN_XLSX_SN_LUT[self._tsv_basename]

    def _get_jsonl_filepath(self) -> str:
        """Get xlsx file sheet name for the TSV file."""
        valid_tsv_suffixes = [".tsv", ".tsv.gz"]
        jsonl_filepath = None
        for val_tsv_sfx in valid_tsv_suffixes:
            if self._tsv_filepath.endswith(val_tsv_sfx):
                jsonl_filepath = "{}.jsonl.gz".format(
                    self._tsv_filepath[:-len(val_tsv_sfx)])

        if jsonl_filepath is None:
            raise ValueError("Invalid TSV file path {}.".format(
                self._tsv_filepath))

        if not os.path.isfile(jsonl_filepath):
            raise ValueError("Corresponding JSONL file {} for TSV file {} "
                             "does not exist.".format(jsonl_filepath,
                                                      self._tsv_filepath))

        return jsonl_filepath

    def _read_tsv(self) -> pd.DataFrame:
        """Read TSV file as pd.DataFrame."""
        # dtype=str reads all columns as strings, because values are not
        # processed in this script
        #
        # na_filter=False reads all values as non-missing, because the missing
        # values are not processed in this script
        tsv_df = pd.read_csv(self._tsv_filepath,
                             sep="\t",
                             dtype=str,
                             na_filter=False)
        return tsv_df

    def _get_col_disp_order(self) -> List[str]:
        """Get column display order"""
        col_order = self._col_names
        if self._tsv_basename == VARIANT_SNV_FN:
            # Move Variant_ID_hg38 and Protein_change columns after the
            # Gene_symbol column
            assert "Variant_ID_hg38" in col_order
            assert "Protein_change" in col_order
            assert "Gene_symbol" in col_order
            left_most_cols = [
                "Gene_symbol", "Variant_ID_hg38", "Protein_change"
            ]
            col_order = (left_most_cols +
                         [x for x in col_order if x not in left_most_cols])
        return col_order

    def _get_col_disp_names(self) -> List[str]:
        """Get column display names"""
        # Change RMTL to PMTL as column display name, according to slack message
        # https://opentargetspediatrics.slack.com/archives/C021ZLY33S7/p1627998314004300
        col_disp_name_lut = {
            "RMTL": "PMTL"
        }
        col_disp_names = [col_disp_name_lut.get(x, x) for x in self._col_names]
        # Replace "_" in column names with " "
        col_disp_names = [x.replace("_", " ") for x in col_disp_names]
        return col_disp_names

    def _get_row_sample_df(self, n_sample_rows: int) -> pd.DataFrame:
        """Get row sample data frame"""
        if n_sample_rows > self._df.shape[0]:
            # If n_sample_rows > number of rows in self._df, take full df and
            # add blank lines
            n_blank_rows = n_sample_rows - self._df.shape[0]
            blank_row_df = pd.DataFrame(index=pd.RangeIndex(start=0,
                                                            stop=n_blank_rows,
                                                            step=1),
                                        columns=self._col_names).fillna("")
            sample_df = pd.concat([self._df.copy(), blank_row_df],
                                  ignore_index=True)
        else:
            # If n_sample_rows <= self._df.shape[0], sample n_sample_rows.
            sample_weights = None
            # Add sample weights to favor rows with non-empty RMTL
            if "RMTL" in self._col_names:
                rmtl_list = self._df["RMTL"].tolist()
                n_empty_rmtl = sum([x == "" for x in rmtl_list])
                n_non_empty_rmtl = sum([x != "" for x in rmtl_list])
                empty_rmtl_weight = 1
                if n_non_empty_rmtl != 0:
                    # no divide by 0
                    non_empty_rmtl_weight = n_empty_rmtl / n_non_empty_rmtl * 2
                else:
                    non_empty_rmtl_weight = 1
                sample_weights = []
                for x in rmtl_list:
                    if x == "":
                        sample_weights.append(empty_rmtl_weight)
                    else:
                        sample_weights.append(non_empty_rmtl_weight)

            sample_df = self._df.sample(n=n_sample_rows,
                                        replace=False,
                                        weights=sample_weights,
                                        random_state=20210811,
                                        axis=0)
            # sort rows by index labels
            sample_df = sample_df.sort_index()

        assert sample_df.shape[0] == n_sample_rows
        return sample_df

    def _get_col_disp_order_name_df(self) -> pd.DataFrame:
        """Get column display order name data frame.

        A column display order name pd.DataFrame has the following rows:

        - Column names in the JSONL/TSV files.
        - Column names for PedOT table view display.
        - 60 sample rows of table values.
        """
        n_sample_rows = 60
        sample_df = self._get_row_sample_df(n_sample_rows)
        # Insert a row at top as display name
        col_disp_name_df = pd.DataFrame([self._get_col_disp_names()],
                                        columns=self._col_names)
        sample_df = pd.concat([col_disp_name_df, sample_df], ignore_index=True)
        # Insert a row at bottom as column name
        # Columns will not be output into xlsx sheet
        col_name_df = pd.DataFrame([self._col_names], columns=self._col_names)
        sample_df = pd.concat([sample_df, col_name_df], ignore_index=True)
        # Reorder columns
        sample_df = sample_df[self._get_col_disp_order()]
        # Insert a column at beginning as row annotation
        sample_df.insert(
            loc=0,
            column="Row annotation",
            value=(
                ["Column display name on PedOT website"] +
                ["Sample row #{}".format(x + 1) for x in range(n_sample_rows)] +
                [("JSON object key name for programming use "
                  "(not displayed on PedOT website)")]))
        # Insert a column at beginning as user guide
        user_guide_column_val = [
            "User guide",
            "The first two columns will not be displayed on PedOT website.",
            ("Change the order of columns to advise PedOT website column "
             "display orders."),
            ("Change the values of the first row to advise PedOT website "
             "column display names."),
            ("Delete one or more columns to advise PedOT website not to "
             "display the deleted columns."),
            ("Please do not change the order or values of the first two "
             "columns."),
            "Please do not delete the first or the second column.",
            "Please do not change the values of the last row.",
            "Please do not delete any row."
        ]  # yapf: disable
        sample_df.insert(
            loc=0,
            column="User guide",
            value=(user_guide_column_val + [
                ""
                for x in range(sample_df.shape[0] - len(user_guide_column_val))
            ]))
        return sample_df

    def write_xlsx_sheet(self, xlsx_writer: pd.ExcelWriter) -> None:
        """Write column display order name data frame into a xlsx sheet"""
        self.col_disp_order_name_df.to_excel(xlsx_writer,
                                             sheet_name=self.xlsx_sheet_name,
                                             header=False,
                                             index=False,
                                             freeze_panes=(1, 2))
