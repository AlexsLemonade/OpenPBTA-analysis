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
                                       fusion_df = NULL,
                                       gene_list = NULL,
                                       color_palette) {
  # Given maf, cnv, and fusion data.frames prepared in `00-map-to-sample_id.R`,
  # along with a list of genes prepared in the `interaction-plots` directory,
  # plot an oncoprint landscape figure.
  #
  # Args:
  #   maf_df: data.frame with data from a MAF file
  #   cnv_df: data.frame with copy number variant data
  #   fusion_df: data.frame with fusion data. This is NULL by default.
  #   gene_list: a list or vector of genes that should be represented on the
  #              oncoprint. This is NULL by default.
  #   color_palette: vector of colors (hex codes) corresponding to the data
  #                  that will be represented on the oncoprint
  
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
  
  if (!is.null(gene_list)) {
    #### Specify genes
    # Get top mutated this data and goi list
    gene_sum <- mafSummary(maf_object)$gene.summary
    
    # Subset for genes in the histology-specific list
    subset_gene_sum <-
      subset(gene_sum, Hugo_Symbol %in% gene_list)
    
    # Get top altered genes
    gene_ordered <-
      subset_gene_sum[order(subset_gene_sum$AlteredSamples, decreasing = TRUE), ]
    
    # Select n top genes
    num_genes <-
      ifelse(nrow(gene_ordered) > 20, 20, nrow(gene_ordered))
    gene_ordered_num <- gene_ordered[1:num_genes, ]
    gene_list <- gene_ordered_num$Hugo_Symbol
  }
  
  #### Plot Oncoprint
  
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
    genes = gene_list,
    logColBar = TRUE,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = TRUE,
    removeNonMutated = TRUE,
    annotationFontSize = 0.7,
    SampleNamefontSize = 0.5,
    fontSize = 0.7,
    colors = color_palette
    ## TODO: Implement use of color palettes in `figures/palettes` directory
  )
  
}