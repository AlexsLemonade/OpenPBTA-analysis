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
                     end_col = "end",
                     y_val = "none",
                     color_col = "none") {
  # Given a data.frame, set it up for use with circlize. Note that start and end
  # columns cannot be the same.
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
  #   color_col: The column name of the colors you would like to use (Optional)
  #
  # Returns:
  #  Data.frame formatted how circlize wants it

  # Filter the samples specified
  if (sample_names[1] != "all") {
    df <- df %>%
      # Filter based on the samples
      dplyr::filter(!!rlang::sym(samples_col) %in% sample_names)
  }

  # Set up the data how circlize wants it
  bed_df <- df %>%
    # Pull out the specified columns
    dplyr::select(
      chr = !!rlang::sym(chr_col),
      start = !!rlang::sym(start_col),
      end = !!rlang::sym(end_col)
    )

  # Don't have to have a y_val
  if (y_val != "none") {
    bed_df <- bed_df %>%
      dplyr::mutate(y_val = as.numeric(dplyr::pull(df, !!rlang::sym(y_val))))
  }
  # Don't have to have a color_col
  if (color_col != "none") {
    bed_df <- bed_df %>%
      dplyr::mutate(color_col = dplyr::pull(df, !!rlang::sym(color_col)))
  }
  # Circlize is particular about what these values need to be
  bed_df <- bed_df %>%
    dplyr::mutate(
      chr = as.character(chr),
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    # Circlize wants things to be only a data.frame NOT a tibble
    as.data.frame()

  # Stop if there isn't data in this data.frame
  if (nrow(bed_df) == 0) {
    stop("There are 0 rows of data. Check your samples_col arg and the data.frame you supplied. ")
  }

  # Check if the chromosomes have the `chr` format. If not, add it.
  if (!grepl("chr", bed_df$chr[1])) {
    bed_df <- bed_df %>%
      dplyr::mutate(chr = paste0("chr", chr))
  }

  # Return formatted data.frame
  return(bed_df)
}


circos_map_plot <- function(df,
                            add_track = FALSE,
                            samples_col = "samples",
                            sample_names = NULL,
                            chr_col = "chrom",
                            start_col = "start",
                            end_col = "end",
                            color_col = "none",
                            y_val = NULL,
                            track_height = .15,
                            type = "point",
                            rect_height = 0.4,
                            cytoband = TRUE,
                            single_color = "none") {
  # Given a data.frame with chromosomal coordinates, and a corresponding data value to plot,
  # make a circos plot or add a circos track to an existing plot.
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
  #   color_col: a column with color specifications for each data point.
  #   y_val: The column name of the value you would like to plot
  #   track_height: a number between 0 and 1 that designates height, 1 = the full diameter of the circos plot.
  #   type: Type of plot the track should be. Options are line, point, rect
  #   rect_height: The added height (plus and minus y_val) that should be plotted.
  #   cytoband: TRUE/FALSE indicating whether you want a cytoband on the outermost
  #             of the plot. Default is TRUE.
  #   single_color: A single color to choose.
  #
  # Returns:
  #  Circos plot (or new circos plot track) of the data provided.

  # Set up the data how circlize wants it
  bed_df <- prep_bed(
    df = df,
    samples_col = samples_col,
    sample_names = sample_names,
    chr_col = chr_col,
    start_col = start_col,
    end_col = end_col,
    y_val = y_val,
    color_col = color_col
  )

  # Establish ymin and ymax before further manipulations
  y_min <- min(bed_df$y_val, na.rm = TRUE)
  y_max <- max(bed_df$y_val, na.rm = TRUE)

  # Can't have identical y_min and y_max, this is just so CircleCI runs even if 
  # the subset data is wonky
  if (y_min == y_max) {
    y_max <- y_max + 0.001 
    warning("ymax and ymin are identical")
  }
    
  # Tell them only one color is allowed
  if (length(single_color) > 1) {
    warning("Only a single color is allowed for the `single_color` argument, 
            only the first item will be used.")
  }

  # If a single color is specified:
  if (single_color != "none" & color_col == "none") {
    # If a single color is specified, replace the color column with that value
    bed_df <- bed_df %>%
      dplyr::mutate(color_col = single_color[1])

    # If both single_color and color_col have been used, spit back an error
  } else if (single_color != "none" & color_col != "none") {
    stop("You've specified a `single_color` and a color_col, need to pick one or 
         the other.")

    # If neither have been used, just use black
  } else if (single_color == "none" & color_col == "none") {
    # If no color key is specified, just make the color black.
    bed_df <- bed_df %>%
      dplyr::mutate(color_col = "black")
  }

  # If we aren't adding on a track to a current plot...
  if (!add_track) {
    # Initialize the plot
    circlize::circos.clear()
    circlize::circos.par(start.degree = 90)

    # Initialize with or without the cytoband
    if (cytoband) {
      circlize::circos.initializeWithIdeogram(species = "hg38")
    } else {
      circlize::circos.initializeWithIdeogram(
        plotType = c("axis", "labels"),
        species = "hg38"
      )
    }
  }

  if (type == "point") {
    # Add scatterplot track
    circlize::circos.genomicTrackPlotRegion(bed_df,
      track.height = track_height,
      # Establish the ylim
      ylim = c(y_min, y_max),
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicPoints(region,
          value[[1]],
          col = value[[2]]
        )
      }
    )
  } else if (type == "line") {
    # Need to separate the lines into their own columns
    spread_df <- bed_df %>%
      # Remove NA values
      dplyr::filter(!is.na(y_val)) %>%
      tibble::rowid_to_column() %>%
      # Make each color its own column
      tidyr::spread(color_col, y_val) %>%
      dplyr::select(-rowid)

    # Palette color
    palette_col <- colnames(spread_df)[4:ncol(spread_df)]

    # Add line track
    circlize::circos.genomicTrackPlotRegion(spread_df,
      track.height = track_height,
      # Establish the ylim
      ylim = c(y_min, y_max),
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicLines(region,
          value,
          col = palette_col
        )
      }
    )
  } else if (type == "rect") {
    # Add a rectangle track
    circlize::circos.genomicTrackPlotRegion(bed_df,
      track.height = track_height,
      # Establish the ylim
      ylim = c(y_min, y_max),
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicRect(region, value[[1]],
          ytop = value[[1]] + rect_height,
          ybottom = value[[1]] - rect_height,
          col = value[[2]],
          border = value[[2]]
        )
      }
    )
  }
}

circos_map_transloc <- function(df,
                                add_track = FALSE,
                                sample_names = NULL,
                                samples_col = "samples",
                                chr_col_1 = "chrom1",
                                chr_col_2 = "chrom2",
                                start_col_1 = "start1",
                                start_col_2 = "start2",
                                end_col_1 = "end1",
                                end_col_2 = "end2",
                                color_col = "none",
                                cytoband = TRUE) {
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
  #   color_col: a column with color specifications for each data point.
  #   cytoband: TRUE/FALSE indicating whether you want a cytoband on the outermost
  #             of the plot. Default is TRUE.
  # Returns:
  #  Circos plot (or new circos plot track) of the translocation data provided.

  # Set up the df_1 how circlize wants it
  bed_df_1 <- prep_bed(
    df = df,
    samples_col = samples_col,
    sample_names = sample_names,
    chr_col = chr_col_1,
    start_col = start_col_1,
    end_col = end_col_1,
    color_col = color_col
  ) # %>%
  # dplyr::mutate(value1 = rnorm(nrow(.)))

  # Set up the df_2 how circlize wants it
  bed_df_2 <- prep_bed(
    df = df,
    samples_col = samples_col,
    sample_names = sample_names,
    chr_col = chr_col_2,
    start_col = start_col_2,
    end_col = end_col_2,
    color_col = color_col
  ) # %>%
  # dplyr::mutate(value1 = rnorm(nrow(.)))

  # If we aren't adding on a track to a current plot...
  if (!add_track) {
    # Initialize the plot
    circlize::circos.clear()
    circlize::circos.par(start.degree = 90)

    # Initialize with or without the cytoband
    if (cytoband) {
      circlize::circos.initializeWithIdeogram(species = "hg38")
    } else {
      circlize::circos.initializeWithIdeogram(
        plotType = c("axis", "labels"),
        species = "hg38"
      )
    }
  }

  if (color_col == "none") {
    # Create random palette if no color palette was called
    color_palette <- circlize::rand_color(nrow(bed_df_1), transparency = 0.5)
  } else {
    # Pull palette if the color_col was called
    color_palette <- dplyr::pull(bed_df_1, !!rlang::sym(color_col))
  }
  circlize::circos.genomicLink(
    bed_df_1,
    bed_df_2,
    col = color_palette,
    border = color_palette
  )
}
