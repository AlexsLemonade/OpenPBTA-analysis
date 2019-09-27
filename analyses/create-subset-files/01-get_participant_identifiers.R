# J. Taroni for CCDL 2019
#
#

# cnv
# fusion
# snv
# sv
# everything that is an RDS file is some kind of expression data
library(readr)
`%>%` <- dplyr::`%>%`

get_biospecimen_ids <- function(filename, id_mapping_df) {
  # Given a supported OpenPBTA file, return the participant IDs corresponding
  # to the biospecimen IDs contained within that file
  #
  # Args:
  #   filename: the full path to a supported OpenPBTA file
  #   id_mapping_df: the data.frame that contains mapping between biospecimen
  #                  IDs and participant IDs
  #
  # Returns:
  #   a vector of unique participant IDs

  message(paste("Reading in", filename, "..."))

  # where the biospecimen IDs come from in each file depends on the file
  # type -- that is why we need all of this logic
  if (grepl("pbta-snv", filename)) {
    # in a column 'Tumor_Sample_Barcode'
    snv_file <- data.table::fread(filename,
                                  skip = 1,  # skip version string
                                  data.table = FALSE)
    biospecimen_ids <- unique(snv_file$Tumor_Sample_Barcode)
  } else if (grepl("pbta-cnv", filename)) {
    # in a column 'ID'
    cnv_file <- read_tsv(filename)
    biospecimen_ids <- unique(cnv_file$ID)
  } else if (grepl("pbta-fusion", filename)) {
    # in a column 'tumor_id'
    fusion_file <- read_tsv(filename)
    biospecimen_ids <- unique(fusion_file$tumor_id)
  } else if (grepl("pbta-sv", filename)) {
    # in a column 'Kids.First.Biospecimen.ID.Tumor'
    sv_file <- data.table::fread(filename, data.table = FALSE)
    biospecimen_ids <- unique(sv_file$Kids.First.Biospecimen.ID.Tumor)
  } else if (grepl(".rds", filename)) {
    # any column name that contains 'BS_' is a biospecimen ID
    expression_file <- read_rds(filename) %>%
      dplyr::select(dplyr::contains("BS_"))
    biospecimen_ids <- unique(colnames(expression_file))
  } else {
    # error-handling
    stop("File type unrecognized by 'get_biospecimen_ids'")
  }

  # map from biospecimen ID to participant IDs and return unique participant IDs
  participant_ids <-  id_mapping_df %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% biospecimen_ids) %>%
    dplyr::pull(Kids_First_Participant_ID)
  return(unique(participant_ids))

}

# take a directory
data_directory <- "data"

# list all files we are interested in subsetting and can support
supported_files_string <-
  "pbta-snv|pbta-cnv|pbta-fusion|pbta-isoform|pbta-sv|pbta-gene"
files_to_subset <- list.files(data_directory,
                              pattern = supported_files_string,
                              full.names = TRUE)
# for testing, get rid of the vardict file
files_to_subset <- files_to_subset[-grepl("vardict", files_to_subset)]

# get the participant ID to biospecimen ID
id_mapping_df <- read_tsv(file.path(data_directory, "pbta-histologies.tsv")) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID) %>%
  dplyr::distinct()

participant_id_list <- purrr::map(files_to_subset,
                                  ~ get_biospecimen_ids(.x, id_mapping_df)) %>%
  stats::setNames(files_to_subset)
