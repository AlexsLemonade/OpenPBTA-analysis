# Set up and calculate functions for handling MAF data
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################

set_up_maf <- function(maf_df, metadata_df = NULL, vaf_cutoff = 0) {
  # Creates these new variables from a MAF formatted data.frame provided: VAF,
  # mutation_id, base_change, change. Optionally can tack on metadata columns
  # which will be matched using a `Tumor_Sample_Barcode` field. Lastly, any
  # columns that contain all `NA` values will be removed.
  #
  # Args:
  #   maf_df: a maf formatted data.frame
  #   metadata: a data.frame with metadata that you would like to merge with the
  #             maf_df and it's newly calculated variables. (Optional)
  #
  # Returns:
  #   a data.frame with all the original information from the MAF data.frame
  #   maf object but with these new variables: VAF, mutation_id, base_change,
  #   change, coding. If metadata_df is specified, then it will also have
  #   those columns added.
  

}

