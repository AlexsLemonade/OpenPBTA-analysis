# Functions for circos plots
#
# C. Savonen for ALSF - CCDL
#
# 2020

circos_map_plot <- function(df = sv_df
                            add_track = FALSE
                            samples_col = "Kids.First.Biospecimen.ID.Tumor"
                            sample_names = "BS_XEVMEYFS"
                            chr_col = "SV.chrom"
                            start_col = "SV.start"
                            end_col =  "SV.end"
                            y_val = "AnnotSV.ranking"
                            track_height = .15
                            type = "point"
                            colour_key =  "SV.type" 
                            palette = "YlGnBl" ) {
  # Given a GRanges object, plot it y value along its chromosomal mappings using
  # ggbio.
  #
  # Args:
  #   df: a data.frame with the chromosomal coordinates and y value to plot
  #   add_track: If true, adds a track to a current plot, if FALSE, starts a new plot.
  #   samples_col: a character string specifying the samples column which can be used to filter by.
  #   sample_names: a character string that specifies values to keep from `samples_col` column. 
  #                 "all" keeps all samples in. 
  #   chr_col: a character string that specifies the chromosomes column name.
  #   start_col: a character string that specifies the start coordinate column name.
  #   end_col: a character string that specifies the end coordinate column name.
  #   y_val: The column name of the value you would like to plot
  #   track_height: a number between 0 and 1 that designates height, 1 = the full diameter of the circos plot. 
  #   type: Type of plot the track should be. Options are line, point, rect
  #   colour_key: a character string designating a column that contains information
  #               to color by. 
  #
  # Returns:
  #  Circos plot of the y value given.
  
  # Filter the samples specified
  if(sample_names != "all") {
    df <- df %>% 
      # Filter based on the samples
      dplyr::filter(!!rlang::sym(samples_col) %in% sample_names) 
  }
  
  # Set up the data how circlize wants it
  bed_df <- df %>% 
    # Pull out the specified columns
    dplyr::select(chr = !!rlang::sym(chr_col), 
                  start = !!rlang::sym(start_col),
                  end = !!rlang::sym(end_col), 
                  value = !!rlang::sym(y_val)) %>% 
    # Circlize is particular about what these values need to be
    dplyr::mutate(chr = as.character(chr), 
                  start = as.integer(start), 
                  end = as.integer(end), 
                  value = as.numeric(value))
  
  # Establish ymin and ymax before further manipulations
  y_min <- min(bed_df$value, na.rm = TRUE)
  y_max <- max(bed_df$value, na.rm = TRUE)
  
  # Stop if there isn't data in this data.frame
  if (nrow(bed_df) == 0){
    stop("There are 0 rows of data. Check your samples_col arg and the data.frame you supplied. ")
  }
  
  # Set up color key if specified
  if (!is.null(colour_key)) {
    bed_df  <- bed_df %>% 
      # Add on this column 
      dplyr::mutate(color_key = as.factor(dplyr::pull(df, !!rlang::sym(colour_key)))) %>% 
      tibble::rowid_to_column() %>% 
      # Make a color column and then spread to each so circlize can use it
      tidyr::spread(color_key, value) %>% 
      dplyr::select(-rowid)
    
    # Set up a palette based on number of factor levels
    palette_col <- topo.colors(ncol(bed_df) - 3)
  }
  
  # Check if the chromosomes have the `chr` format. If not, add it. 
  if(!grepl("chr", bed_df$chr[1])) {
    bed_df <- bed_df %>% 
      dplyr::mutate(chr = paste0("chr", chr))
    }
  
  if (!add_track) {
  # Initialize the plot
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90)
  circlize::circos.initializeWithIdeogram(species = "hg38")
  }
  
  if (type == "point") {
  # Add scatterplot track
  circlize::circos.genomicTrackPlotRegion(bed_df, 
                                          track.height = track_height,
                                          # Establish the ylim
                                          ylim = c(y_min, y_max),
                                          panel.fun = function(region, value, ...) {
                                              circlize::circos.genomicPoints(region, value, 
                                                                             col = palette_col)
                                          })
  } else if (type == "line") {
  # Add line track
  circlize::circos.genomicTrackPlotRegion(bed_df, 
                                          track.height = track_height,
                                          # Establish the ylim
                                          ylim = c(y_min, y_max),
                                          panel.fun = function(region, value, ...) {
                                            circlize::circos.genomicLines(region, value, 
                                                                           col = palette_col)
                                          })
  } else if (type == "rect") {
    # Add the actual track of scatterplots
    circlize::circos.genomicTrackPlotRegion(bed_df, 
                                            track.height = track_height,
                                            # Establish the ylim
                                            ylim = c(y_min, y_max),
                                            panel.fun = function(region, value, ...) {
                                              circlize::circos.genomicRect(region, value, 
                                                                           col = palette_col)
                                            })
  }
}

