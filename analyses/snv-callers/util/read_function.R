# Function for reading TSV or RDS
#
# C. Savonen for ALSF - CCDL
# 2019
#
################################################################################
read_tsv_or_rds <- function(file_name) {
  # Reading a files whether they are TSV or RDS, give an error if it is neither.
  #
  # Args:
  #   file_name: a file you suspect is either a TSV or RDS file that you want
  #              to read in.

  # Return: read in file

  # What kind of file is this?
  file_suffix <- substr(file_name, nchar(file_name) - 2, nchar(file_name))

  # Make to lower in case upper RDS was used
  file_suffix <- tolower(file_suffix)

  # Read in if its TSV
  if (file_suffix  == "tsv") {
    return(readr::read_tsv(file))
  }
  # Read in if its RDS
  if (file_suffix == "rds") {
    return(readr::read_rds(file_name))
  } else {
    # If neither, stop
    stop(paste(file_name, "cannot be read", file_suffix, "is not supported"))
  }
}
