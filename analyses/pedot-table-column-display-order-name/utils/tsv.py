import pandas as pd
import os


class TSV:
    """TSV table for generating PedOT column display order and name xlsx sheet.

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
        self._read_tsv()
        self._set_xlsx_sheet_name()
        self._set_col_disp_order_name_df()

    def _set_xlsx_sheet_name(self) -> None:
        """Set xlsx file sheet name for the TSV file."""
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

        self._jsonl_filepath = jsonl_filepath
        self.xlsx_sheet_name = os.path.basename(jsonl_filepath)

    def _read_tsv(self) -> None:
        """Read TSV file as pd.DataFrame."""
        # dtype=str reads all columns as strings, because values are not
        # processed in this script
        #
        # na_filter=False reads all values as non-missing, because the missing
        # values are not processed in this script
        self._df = pd.read_csv(self._tsv_filepath,
                               sep="\t",
                               dtype=str,
                               na_filter=False)

    def _set_col_disp_order_name_df(self) -> None:
        """Set column display order name data frame.

        A column display order name pd.DataFrame has the following rows:

        - Column names in the JSONL/TSV files.
        - Column names for PedOT table view display.
        - 10 sample rows of table values.
        """
        n_sample_rows = 10
        sample_df = self._df.sample(n=n_sample_rows,
                                    replace=False,
                                    random_state=20210811,
                                    axis=0)
        # Insert a row at top as display name
        col_names = sample_df.columns.tolist()
        col_disp_names = [x.replace("_", " ") for x in col_names]
        col_disp_name_df = pd.DataFrame([col_disp_names], columns=col_names)
        sample_df = pd.concat([col_disp_name_df, sample_df], ignore_index=True)
        # Insert a row at bottom as column name
        # Columns will not be output into xlsx sheet
        col_name_df = pd.DataFrame([col_names], columns=col_names)
        sample_df = pd.concat([sample_df, col_name_df], ignore_index=True)
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
            ("Please do not change the order or values of the first two "
             "columns."),
            "Please do not change the values of the last row."
        ] # yapf: disable
        sample_df.insert(
            loc=0,
            column="User guide",
            value=(user_guide_column_val + [
                ""
                for x in range(sample_df.shape[0] - len(user_guide_column_val))
            ]))
        self.col_disp_order_name_df = sample_df

    def write_xlsx_sheet(self, xlsx_writer: pd.ExcelWriter) -> None:
        self.col_disp_order_name_df.to_excel(xlsx_writer,
                                             sheet_name=self.xlsx_sheet_name,
                                             header=False,
                                             index=False,
                                             freeze_panes=(1, 2))
