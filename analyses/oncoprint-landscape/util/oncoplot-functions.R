# This script defines custom functions to be sourced in the
# `01-plot-oncoprint.R` script of this module.
#
# Chante Bethell for CCDL 2020
#
# # #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("util", "oncoplot-functions.R"))

prepare_maf_object <- function(maf_df,
                               cnv_df,
                               metadata,
                               fusion_df = NULL) {
  # Given maf, cnv, and fusion data.frames prepared in `00-map-to-sample_id.R`,
  # along with the relevant metadata prepared in the `interaction-plots`
  # directory, return a maf object. Note that the fusion data.frame argument is
  # optional here and is NULL by default.
  #
  # Args:
  #   maf_df: data.frame with data from a MAF file
  #   cnv_df: data.frame with copy number variant data
  #   metadata: data.frame with the relavant metadata
  #   fusion_df: data.frame with fusion data. This is NULL by default.
  
  if (!is.null(fusion_df)) {
    # Bind rows of maf and fusion data frames
    maf_df <- dplyr::bind_rows(maf_df, fusion_df)
  }
  
  # Convert into MAF object
  maf_object <-
    read.maf(
      maf = maf_df,
      clinicalData = metadata,
      cnTable = cnv_df,
      removeDuplicatedVariants = FALSE,
      vc_nonSyn = c(
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Splice_Site",
        "Nonsense_Mutation",
        "Nonstop_Mutation",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Missense_Mutation",
        "Fusion",
        "Multi_Hit",
        "Multi_Hit_Fusion",
        "Hom_Deletion",
        "Hem_Deletion",
        "amplification",
        "gain",
        "loss"
      )
    )
  
  # Return the maf object
  return(maf_object)
  
}
