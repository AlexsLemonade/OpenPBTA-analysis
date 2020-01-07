# Functions for chromosomal instability calculations
#
# C. Savonen for ALSF - CCDL
#
# 2020

make_granges <- function(break_df = cnv_breaks, 
                         sample_id = NULL, 
                         samples_col = "samples", 
                         chrom_col =  "chrom", 
                         start_col = "start", 
                         end_col = "end") {
  # For a given breaks data.frame make a GenomicRanger object from it. Optionally
  # can filter to a single samples' data. 
  #
  # Args:
  #   break_df: for a data.frame with chromosomal coordinates and sample IDs, 
  #             any other columns in this data.frame will also be carried along. 
  #   sample_id: sample ID for which the data needs to be extracted and made into 
  #              a GenomicRanges object. If "all" is designated, all samples will
  #              be kept. 
  #   samples_col: character string that indicates the column name with the 
  #                sample ID information. Default is "samples".
  #   chrom_col: character string that indicates the column name with the 
  #              chromosome information. Default is "chrom".
  #   start_col: character string that indicates the column name with the 
  #              start coordinate. Default is "start".
  #   end_col: character string that indicates the column name with the 
  #            end coordinate. Default is "end".
  #
  # Returns:
  # A Genomic Ranges formatted object for the particular sample that has the 
  # break points from the data.frame
  
  # List the columns we need
  needed_cols <- c(samples_col, chrom_col, start_col, end_col)
  
  # Get logical vector indicating which are in metadata
  found_cols <- (needed_cols %in% colnames(break_df))
  
  # If not all the columns are found, stop
  if (!all(found_cols)) {
    stop(cat(
      "The following column names specified for making the GenomicRanges object: \n",
      paste(needed_cols[which(!found_cols)], collapse = "\n"),
      "\n ...were not found in the specified metadata data.frame.",
      "Check your `_col` arguments."
    ))
  }
  
  if (sample_id != "all") {
  # Extract samples' info from the data.frame
    break_df <- break_df %>% 
      dplyr::filter(!!rlang::sym(samples_col) == sample_id)
    
    # Stop if the sample doesn't exist in the data.frame
    if (nrow(break_df) == 0) {
      stop("No sample data was found for the given `sample_id` argument.")
    }
  }

  # Make GRanges for CNV data
  granges <- GenomicRanges::GRanges(
    seqnames = dplyr::pull(break_df, !!rlang::sym(chrom_col)),
    ranges = IRanges::IRanges(
      start = dplyr::pull(break_df, !!rlang::sym(start_col)),
      end = dplyr::pull(break_df, !!rlang::sym(end_col))
    ),
    mcols = break_df 
  )
  
  return(granges)
}

break_density <- function(sv_breaks = NULL, 
                          cnv_breaks = NULL, 
                          sample_id = NULL,
                          window_size = 1e6,
                          max_gap = 0,
                          chr_sizes_list = NULL,
                          samples_col_cnv = "samples", 
                          chrom_col_cnv =  "chrom", 
                          start_col_cnv = "start", 
                          end_col_cnv = "end", 
                          samples_col_sv = "samples", 
                          chrom_col_sv =  "chrom", 
                          start_col_sv = "start", 
                          end_col_sv = "end"
                          ) {
  # For a given breaks data.frame calculate the breaks density for a sliding 
  # window across the genome. Returns the data as a GenomicRanges object for 
  # easy mapping with ggbio. 
  #
  # Args:
  #   sv_breaks: a data.frame with the breaks for the SV data. 
  #   cnv_breaks: a data.frame with the breaks for the SV data.
  #   sample_id: sample ID for which the data needs to be extracted and made into 
  #              a GenomicRanges object. If "all" is designated, all samples will
  #              be kept. 
  #   window_size: What size windows to calculate break density. Default is 1kb. 
  #   samples_col_sv/cnv: character string that indicates the column name with the 
  #                sample ID information. Default is "samples". Will be passed to
  #                make_granges function. 
  #   chrom_col_sv/cnv: character string that indicates the column name with the 
  #              chromosome information. Default is "chrom". Will be passed to
  #              make_granges function. 
  #   start_col_sv/cnv: character string that indicates the column name with the 
  #              start coordinate. Default is "start". Will be passed to
  #              make_granges function. 
  #   end_col_sv/cnv: character string that indicates the column name with the 
  #            end coordinate. Default is "end". Will be passed to
  #            make_granges function. 
  #
  # Check that a sample ID has been specified. 
  if (is.null(sample_id)) {
    stop("No sample ID has been specified. Use the `sample_id` argument.")
  }
  # If no data.frame is provided, stop
  if (is.null(cnv_breaks) && is.null(sv_breaks)) {
    stop("Need at least CNV or SV breaks data.frame provided.")
  }
  # Make CNV into GRanges
  if (!is.null(cnv_breaks)) {
    cnv_ranges <- make_granges(cnv_breaks, 
                              sample_id = sample_id, 
                              samples_col = samples_col_cnv,
                              chrom_col =  chrom_col_cnv, 
                              start_col = start_col_cnv, 
                              end_col = end_col_cnv) 
    # This will get written over if there are both datasets
    combo_ranges <- cnv_ranges
  }
  # Make SV into GRanges
  if (!is.null(sv_breaks)) {
    sv_ranges <- make_granges(sv_breaks, 
                             sample_id = sample_id,
                             samples_col = samples_col_sv, 
                             chrom_col =  chrom_col_sv, 
                             start_col = start_col_sv, 
                             end_col = end_col_sv) 
    # This will get written over if there are both datasets
    combo_ranges <- sv_ranges
  }
  
  # If both datasets are given, reduce them into one. 
  if (!is.null(cnv_breaks) && !is.null(sv_breaks)) {
    combo_ranges <- GenomicRanges::union(cnv_ranges, 
                                         sv_ranges)
  } 
  
  # Create genome bins
  bins <- GenomicRanges::tileGenome(chr_sizes_list, 
                                    tilewidth = window_size)

  # Uncompress GRangesList
  bins <- unlist(bins)
  
  # Get counts for each genome bin
  bin_counts <- GenomicRanges::countOverlaps(bins, 
                                             combo_ranges, 
                                             max_gap)
  
  # Store count info
  bins@elementMetadata@listData$counts <- bin_counts
  
  # Calculate and store density
  bins@elementMetadata@listData$density <- bin_counts/window_size
  
  # Return the GRanges object for mapping purposes
  return(bins)
}

map_density_plot <- function(granges, 
                             y_val, 
                             y_lab, 
                             color,
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
  #  
  # For setting the scale later, need to get y's max
  max_y <- max(
    dplyr::pull(data.frame(granges@elementMetadata@listData), 
                  !!rlang::sym(y_val))
    )
  # Make the density plot
  density_plot <- ggbio::autoplot(granges, ggplot2::aes(y = !!rlang::sym(y_val)),
                                  geom = "line", scales = "free_x", space = "free_x", 
                                  colour = color) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 3, angle = 45, hjust = 1)) +
    ggplot2::ylab(y_lab) + 
    ggplot2::ggtitle(main_title) +
    ggplot2::scale_y_continuous(breaks = seq(0, max_y, by = 2))
                                
  # Print out plot
  density_plot@ggplot
}

