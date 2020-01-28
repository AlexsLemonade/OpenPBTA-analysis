# Functions for circos plots
#
# C. Savonen for ALSF - CCDL
#
# 2020
# 
prep_bed <- function(df,
                     samples_col = "samples",
                     sample_names = "all",
                     chr_col = "chrom",
                     start_col = "start",
                     end_col =  "end",
                     y_val = "none") {  
  # Given a data.frame, set it up for use with circlize. Note that start and end
  #  columns cannot be the same. 
  #
  # Args:
  #   df: a data.frame with the chromosomal coordinates and y value to plot
  #   samples_col: a character string specifying the samples column which can be used to filter by.
  #   sample_names: a character string that specifies values to keep from `samples_col` column. 
  #                 "all" keeps all samples in. 
  #   chr_col: a character string that specifies the chromosomes column name.
  #   start_col: a character string that specifies the start coordinate column name. 
  #   end_col: a character string that specifies the end coordinate column name.
  #   y_val: The column name of the value you would like to plot (Optional)
  #
  # Returns:
  #  Data.frame formatted how circlize wants it 
  
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
                  end = !!rlang::sym(end_col))

  # Don't have to have a y_val
  if (y_val != "none") {
    bed_df <- bed_df %>% 
      dplyr::mutate(y_val = as.numeric(dplyr::pull(df, !!rlang::sym(y_val))))
  } 
  # Circlize is particular about what these values need to be
  bed_df <- bed_df %>% 
    dplyr::mutate(chr = as.character(chr), 
                start = as.integer(start), 
                end = as.integer(end)) %>% 
    # Circlize wants things to be only a data.frame NOT a tibble
    as.data.frame()
  
  # Stop if there isn't data in this data.frame
  if (nrow(bed_df) == 0){
    stop("There are 0 rows of data. Check your samples_col arg and the data.frame you supplied. ")
  }
  
  # Check if the chromosomes have the `chr` format. If not, add it. 
  if(!grepl("chr", bed_df$chr[1])) {
    bed_df <- bed_df %>% 
      dplyr::mutate(chr = paste0("chr", chr))
  }
  
  # Return formatted data.frame
  return(bed_df)
}

df = cnv_df 
add_track = FALSE
samples_col = "ID"
sample_names = "BS_007JTNB8"
chr_col = "chrom"
start_col = "loc.start"
end_col =  "loc.end"
y_val = "copy.num"
track_height = .15
type = "rect"
colour_key = "direction"
palette = "YlGnBu"

circos_map_plot <- function(df,
                            add_track = FALSE,
                            samples_col = "samples",
                            sample_names = NULL,
                            chr_col = "chrom",
                            start_col = "start",
                            end_col =  "end",
                            y_val = NULL,
                            track_height = .15,
                            type = "rect",
                            colour_key = NULL, 
                            palette = "YlGnBu") {
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
  #               to color by. Can be numeric or a factor.  
  #   palette: the color brewer palette you would like to be used. Run `RColorBrewer::display.brewer.all()` to see options. 
  #
  # Returns:
  #  Circos plot (or new circos plot track) of the data provided. 
  
  # Filter the samples specified
  if (sample_names != "all") {
    df <- df %>% 
      # Filter based on the samples
      dplyr::filter(!!rlang::sym(samples_col) %in% sample_names) 
  }
  
  # Set up the data how circlize wants it
  bed_df <- prep_bed(df = df,
                     samples_col = samples_col,
                     sample_names = sample_names,
                     chr_col = chr_col,
                     start_col = start_col,
                     end_col = end_col,
                     y_val = y_val)
  
  # Establish ymin and ymax before further manipulations
  y_min <- min(bed_df$y_val, na.rm = TRUE)
  y_max <- max(bed_df$y_val, na.rm = TRUE)

  # Set up color key if specified
  if (!is.null(colour_key)) {
    # Add on color_key as a new column
    bed_df  <- bed_df %>% 
      dplyr::mutate(color_key = dplyr::pull(df, !!rlang::sym(colour_key))) %>%
      tibble::rowid_to_column() 
    
    # If color key is a factor...
    if (is.factor(bed_df$color_key)) {
      # Determine the number of levels
      n_levels <- length(levels(bed_df$color_key))
      
      # ColorBrewer complains if given a number less than 3 
      # Set up a palette based on number of factor levels
      palette_col <- RColorBrewer::brewer.pal(ifelse(n_levels < 3, 3, n_levels), 
                                              name = palette)
      
      # Make a recode color key list 
      level_key <- palette_col[1:n_levels]
      
      # Make the palette colors the names in this list
      names(level_key) <- levels(bed_df$color_key)
      
      # Set up the color ready data.frame
      bed_df  <- bed_df %>% 
        # Make column that specifies the color for each value
        dplyr::mutate(color_key = dplyr::recode(color_key, !!!level_key), 
                      color_key = as.character(color_key))
      
    # If color key is numeric...
    } else if (is.numeric(bed_df$color_key)) {
      # Make color palette based on 5 colors
      palette_col <- RColorBrewer::brewer.pal(5, name = palette)
      
      # Make color ramp function based on quantiles and specified palette
      color_fun <- circlize::colorRamp2(breaks = quantile(bed_df$color_key, 
                                                            c(0.15, 0.35, 0.5, 0.65, 0.85)), 
                                        colors = palette_col)
      # Set up the color ready data.frame
      bed_df  <- bed_df %>% 
        # Make column that specifies the color for each value
        dplyr::mutate(color_key = color_fun(value))
    } else {
      stop("`colour_key` argument can only be a numeric or factor variable. It is neither. 
           Check the column data you specified for `colour_key`.")
    }
  } else {
    # If no color key is specified, just make the color black.
    palette_col <- "black"
  }
  
  # Check if the chromosomes have the `chr` format. If not, add it. 
  if(!grepl("chr", bed_df$chr[1])) {
    bed_df <- bed_df %>% 
      dplyr::mutate(chr = paste0("chr", chr))
    }
  
  # If we aren't adding on a track to a current plot...
  if (!add_track) {
  # Initialize the plot
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90)
  circlize::circos.initializeWithIdeogram(species = "hg38")
  }
  
  if (type == "point") {
    # Circlize wants things in their own column for each color
    bed_df  <- bed_df %>% 
      # Make a color column and then spread to each so circlize can use it
      tidyr::spread(color_key, value) %>% 
      dplyr::select(-rowid)
    
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
    # For line and point, circlize wants things in their own column for each color
    bed_df  <- bed_df %>% 
      # Make a color column and then spread to each so circlize can use it
      tidyr::spread(color_key, value) %>% 
      dplyr::select(-rowid)
    
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
    # Get rid of rowid
    bed_df <- bed_df %>% 
      dplyr::select(-rowid)
    # Add a rectangle track
    circlize::circos.genomicTrackPlotRegion(bed_df, 
                                            track.height = track_height,
                                            # Establish the ylim
                                            ylim = c(y_min, y_max),
                                            panel.fun = function(region, value, ...) {
                                              i = circlize::getI(...)
                                                circlize::circos.genomicRect(region, value[[1]],
                                                                             ytop = value[[1]] + 0.4, 
                                                                             ybottom = value[[1]] - 0.4, 
                                                                             col = value[[2]], 
                                                                             border = value[[2]])
                                            })
  }
}

df = transloc_df 
add_track = FALSE 
sample_names = "BS_007JTNB8"
samples_col = "biospecimen_id1"
chr_col_1 = "chrom1"
chr_col_2 = "chrom2"
start_col_1 = "start1"
start_col_2 = "start2"
end_col_1 =  "end1"
end_col_2 =  "end2"
y_val_2 = 

circos_map_transloc <- function(df, 
                                add_track = FALSE,
                                sample_names = NULL,
                                samples_col = "samples",
                                chr_col_1 = "chrom1",
                                chr_col_2 = "chrom2",
                                start_col_1 = "start1",
                                start_col_2 = "start2",
                                end_col_1 =  "end1",
                                end_col_2 =  "end2") {
  # Given two BED data.frames, plot translocations on a new or already exisiting circos plot. 
  #
  # Args:
  #   df: a data.frame with two sets of chromosomal coordinates indicating the translocations
  #   add_track: If true, adds a track to a current plot in the global environment, if FALSE, starts a new plot.
  #   sample_names: a character string that specifies values to keep from `samples_col` column. 
  #                 "all" keeps all samples in. 
  #   samples_col: a character string specifying the samples column which can be used to filter by.
  #   chr_col_1/2: a character string that specifies the chromosomes column name for the from and to coordinates respectively.
  #   start_col_1/2: a character string that specifies the start coordinate column name for the from and to coordinates respectively.
  #   end_col_1/2: a character string that specifies the end coordinate column name for the from and to coordinates respectively.
  #   y_val_1/2: The column name of the value you would like to plot for the from and to coordinates respectively.
  #   palette: the color brewer palette you would like to be used. Run `RColorBrewer::display.brewer.all()` to see options. 
  #
  # Returns:
  #  Circos plot (or new circos plot track) of the data provided. 

  # Set up the df_1 how circlize wants it
  bed_df_1 <- prep_bed(df = df,
                       samples_col = samples_col,
                       sample_names = sample_names,
                       chr_col = chr_col_1,
                       start_col = start_col_1,
                       end_col = end_col_1)
  
  # Set up the df_2 how circlize wants it
  bed_df_2 <- prep_bed(df = df,
                      samples_col = samples_col,
                      sample_names = sample_names,
                      chr_col = chr_col_2,
                      start_col = start_col_2,
                      end_col = end_col_2) 
  
  # If we aren't adding on a track to a current plot...
  if (!add_track) {
    # Initialize the plot
    circlize::circos.clear()
    circlize::circos.par(start.degree = 90)
    circlize::circos.initializeWithIdeogram(species = "hg38")
  }
  
  circlize::circos.genomicLink(bed_df_1,
                               bed_df_2, 
                               col = circlize::rand_color(nrow(bed_df_1), transparency = 0.5), 
                               border = NA)
}
