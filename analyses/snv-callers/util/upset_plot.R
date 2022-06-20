# Custom function for upset plots
#
# 2020
# CCDL - C. Savonen

upset_png <- function(detect_mat, plot_file_path, subset_vector = NULL, has_vardict = TRUE) {
  # Given a detection data frame, where each column represents a caller and contains a TRUE/FALSE 
  # for whether a mutation was detected, create an upset plot.
  #
  # Arguments: 
  # detect_mat: a matrix, where each column represents a caller and contains a TRUE/FALSE 
  #            for whether a mutation was detected. 
  # plot_file_path: File path to where the plot will be saved to as png. 
  # subset_vector: (Optional) a vector to be supplied to dplyr::filter() that will subset the
  #                detect_mat and create the upset plot from that subset.
  # has_vardict: Does detect_mat contain VarDict as well? TRUE/FALSE. Default is TRUE.
  #
  # Output: An upset plot that represents the number of mutations detected by the combinations of callers
  
  if (!is.null(subset_vector[1])) {
    # If subset index is set
    detect_mat <- detect_mat[subset_vector, ]
  }
  
  # Set up a list how UpSetR wants it
  detect_list <- list(
    lancet = which(detect_mat[, "VAF_lancet"]),
    mutect = which(detect_mat[, "VAF_mutect"]),
    strelka = which(detect_mat[, "VAF_strelka"])
  )
  
  # Add Vardict if its there
  if (has_vardict) {
    detect_list[["vardict"]] <- which(detect_mat[, "VAF_vardict"])
  }
  
  # Save to PNG
  png(file.path(plot_file_path), width = 1300, height = 900);
  print(
    UpSetR::upset(
      UpSetR::fromList(detect_list), 
      order.by = "freq",
      text.scale = 2,
      point.size = 4,
      mainbar.y.label = "")
  );
  dev.off()
}
