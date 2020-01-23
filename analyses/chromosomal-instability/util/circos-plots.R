# Functions for circos plots
#
# C. Savonen for ALSF - CCDL
#
# 2020

sample_breakpoints <- readr::read_rds(
  file.path(scratch_dir, "sample_breakpoint_densities.RDS"))

circos_map_plot <- function(data = cnv_df,
                            samples_col = "ID",
                            chr_col = "chrom",
                            start_col = "loc.start",
                            end_col =  "loc.end",
                            y_val = "seg.mean",
                            main_title) {
  # Given a GRanges object, plot it y value along its chromosomal mappings using
  # ggbio.
  #
  # Args:
  #   granges: A Granges object to plot
  #   y_val: a character string of the columnname in listData spot of the
  #          GenomicRanges to plot on the y axis
  #   color: a color parameter
  #   y_lab: a character string to use for the ylabel. Will be passed to
  #          ggplot2::ylab argument.
  #   main_title: a character string to use for the main title.
  #
  # Returns:
  #  ggplot of chromosomal mapping of the y value given.
  
  bed_df <- data %>% 
    # Pull out the specified columns
    dplyr::select(chr = !!rlang::sym(chr_col), 
                  start = !!rlang::sym(start_col),
                  end = !!rlang::sym(end_col), 
                  value1 = !!rlang::sym(y_val)) %>% 
    # Circlize is particular about what these values need to be
    dplyr::mutate(chr = as.character(chr), 
                  start = as.integer(start), 
                  end = as.integer(end), 
                  value1 = as.numeric(value1))
  
  # Initialize the plot
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90)
  circlize::circos.initializeWithIdeogram(species = "hg38")
  
  # Add the actual track of scatterplots
  circlize::circos.genomicTrackPlotRegion(bed_df, ylim = c(min(bed_df$value1), max(bed_df$value1)),
                                          panel.fun = function(region, value, ...) {
                                            circlize::circos.genomicPoints(region, value, color = "black")
                                          })
  circlize::circos.clear()
}

