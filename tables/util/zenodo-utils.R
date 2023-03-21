# S. Spielman for CCDL, 2023
# This script holds a utility function for writing figure data CSV tables
#  to the `tables/zenodo-upload/` (or other) directory.


# Define output directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
zenodo_output_dir <- file.path(root_dir, )

#' Write out a CSV-formatted table, ordered by sample, for Zenodo upload. Nothing is returned.
#'
#' @param df Data frame to write to file
#' @param sample_id_column Column containing sample id information to arrange on
#' @param output_filename Output file name. The extension should be `.csv` or
#'   `.csv.gz` if compression is needed
#' @param output_dir Output directory to write to. By default, this is `tables/zenodo-upload/`.
write_zenodo_table <- function(df,
                               sample_id_column,
                               output_filename,
                               output_dir = zenodo_output_dir) {
  # Ensure output_dir exists
  if (!(dir.exists(output_dir))) {
    stop("Error: Output directory provided to `write_zenodo_table()` does not exist.")
  }

  # Order data frame by sample information
  df <- dplyr::arrange(df, {{sample_id_column}})

  # Form final path
  output_filename <- file.path(output_dir, output_filename)

  # Write
  readr::write_csv(df, output_filename)
}