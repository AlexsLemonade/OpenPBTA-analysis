# This script defines custom functions to be sourced in the
# `01-plot-oncoprint.R` script of this module.
#
# Chante Bethell for CCDL 2020
#
# # #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("util", "oncoplot-functions.R"))

prepare_and_plot_oncoprint <- function(maf_df,
                                       cnv_df,
                                       metadata,
                                       fusion_df = NULL,
                                       goi_list = NULL,
                                       color_palette,
                                       running_in_figures_script = FALSE) {
  # Given maf, cnv, and fusion data.frames prepared in `00-map-to-sample_id.R`,
  # along with a list of genes prepared in the `interaction-plots` directory,
  # plot an oncoprint landscape figure.
  #
  # Args:
  #   maf_df: data.frame with data from a MAF file
  #   cnv_df: data.frame with copy number variant data
  #   metadata: data.frame with the relavant metadata
  #   fusion_df: data.frame with fusion data. This is NULL by default.
  #   goi_list: a list or vector of genes of interest that should be
  #             represented on the oncoprint. This is NULL by default.
  #   color_palette: vector of colors (hex codes) corresponding to the data
  #                  that will be represented on the oncoprint
  #   running_in_figures_script: TRUE or FALSE statement indicating whether
  #                              or not this function is being run in the
  #                              figures generation script
  
  if (!is.null(fusion_df)) {
    # Bind rows of maf and fusion data frames
    maf_df <- dplyr::bind_rows(maf_df, fusion_df)
  }
  
  if (!is.null(goi_list)) {
    # Filter `goi_list` to include only the unique genes of interest
    # that are also in `maf_df`
    goi_list <- unique(goi_list[goi_list %in% maf_df$Hugo_Symbol])
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
  
  # If this function is being run in the figures generation script, return the
  # maf object (instead of the final plot)
  if (running_in_figures_script) {
    return(maf_object)
  }
  
  #### Plot Oncoprint
  if (!(running_in_figures_script)) {
    # Given a maf file, plot an oncoprint of the variants in the
    # dataset and save as a png file.
    oncoplot(
      maf_object,
      clinicalFeatures = c(
        "broad_histology",
        "short_histology",
        "reported_gender",
        "tumor_descriptor",
        "molecular_subtype"
      ),
      genes = goi_list,
      logColBar = TRUE,
      sortByAnnotation = TRUE,
      showTumorSampleBarcodes = TRUE,
      removeNonMutated = TRUE,
      annotationFontSize = 0.7,
      SampleNamefontSize = 0.5,
      fontSize = 0.7,
      colors = color_palette
    )
  }
}