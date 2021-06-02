# Function for splitting up MNVs to SNVs
#
# J. Shapiro for ALSF - CCDL
# 2019
#
############################## Custom Function #################################
#' Split multinucleotide variants into single nucleotide calls
#'
#' @param mnv_tbl a table containing MNVs (may be from an sql connection)
#'
#' @return a data frame of SNVs derived from the MNVs
#'   (does not contain any of the original SNVs)
split_mnv <- function(mnv_tbl) {
  mnv_df <- mnv_tbl %>%
    dplyr::filter(Variant_Type %in% c("DNP", "TNP", "ONP")) %>%
    as.data.frame() %>%
    # add a mnv_id for calculating positions, and for potential later
    # reconstitution of MNVs.
    dplyr::mutate(mnv_id = dplyr::row_number()) %>%
    # lancet adds a base to the start of MNV `Allele` fields, so check for that
    # and remove the extra base.
    dplyr::mutate(
      Allele = dplyr::case_when(
        stringr::str_length(Allele) == stringr::str_length(Reference_Allele) ~
        Allele,
        stringr::str_length(Allele) == stringr::str_length(Reference_Allele) + 1 ~
        stringr::str_sub(Allele, 2)
      )
    )
  # separate multibase calls into individual rows
  # Some tables have the `Norm_Seq_Allele`s, so we will preserve and split those
  # if necessary
  if ("Match_Norm_Seq_Allele1" %in% colnames(mnv_df)) {
    mnv_df <- mnv_df %>%
      tidyr::separate_rows(
        Reference_Allele,
        Tumor_Seq_Allele1,
        Tumor_Seq_Allele2,
        Match_Norm_Seq_Allele1,
        Match_Norm_Seq_Allele2,
        Allele,
        sep = "(?<=[A-Za-z])"
      )
  } else {
    mnv_df <- mnv_df %>%
      tidyr::separate_rows(
        Reference_Allele,
        Tumor_Seq_Allele1,
        Tumor_Seq_Allele2,
        Allele,
        sep = "(?<=[A-Za-z])"
      )
  }
  mnv_df <- mnv_df %>%
    # character separation leaves an extra blank
    dplyr::filter(
      Reference_Allele != "",
      Allele != ""
    ) %>%
    dplyr::group_by(mnv_id) %>%
    dplyr::mutate(
      mnv_pos = dplyr::row_number(),
      Start_Position = Start_Position + mnv_pos - 1,
      End_Position = Start_Position
    ) %>%
    dplyr::ungroup()
  return(mnv_df)
}
