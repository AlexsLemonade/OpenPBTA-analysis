# S. Spielman for CCDL, 2023
# This script holds a utility function for writing figure data CSV tables


#' Write out a CSV-formatted table, ordered by sample, for Zenodo upload. Nothing is returned.
#'
#' @param df Data frame to write to file
#' @param sample_id_column Column containing sample id information to arrange on
#' @param output_filepath Output file path. The file extension should be `.csv` or
#'   `.csv.gz` if compression is needed
#' @param output_dir Output directory to write to. By default, this is `tables/zenodo-upload/`.
write_zenodo_table <- function(df,
                               sample_id_column,
                               output_filepath) {

  # Ensure directory exists
  if (!(dir.exists(dirname(output_dir)))) {
    stop("Error: Output directory provided to `write_zenodo_table()` does not exist.")
  }

  # Ensure proper extension
  if (!(stringr::str_ends(output_filepath, "\\.csv|\\.csv\\.gz"))) {
    stop("Error: File name must end in `.csv` or .`csv.gz` for compression.")
  }

  # Order data frame by sample information
  df <- dplyr::arrange(df, {{sample_id_column}})

  # Write
  readr::write_csv(df, output_filename)
}
