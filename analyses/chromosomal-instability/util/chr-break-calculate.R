# Functions for chromosomal instability plots
#
# C. Savonen for ALSF - CCDL
#
# 2020

make_granges <- function(break_df = NULL,
                         sample_id = NULL,
                         samples_col = "samples",
                         chrom_col = "chrom",
                         start_col = "start",
                         end_col = "end") {
  # For a given breaks data.frame make a GenomicRanges object from it. Can
  # filter to a single samples' data.
  #
  # Args:
  #   break_df: for a data.frame with chromosomal coordinates and sample IDs.
  #   sample_id: a character string that designates which data needs to be
  #              extracted and made into a GenomicRanges object by matching the
  #              information in the designated sample_col. If "all" is designated,
  #              all samples will be kept. Multiple samples can be designated
  #              as a character vector.
  #   samples_col: character string that indicates the column name with the
  #                sample ID information. Default is "samples". This information
  #                will be store in the GenomicRanges object in
  #                `elementData@listData$samples`
  #   chrom_col: character string that indicates the column name with the
  #              chromosome information. Default is "chrom".
  #   start_col: character string that indicates the column name with the
  #              start coordinate. Default is "start".
  #   end_col: character string that indicates the column name with the
  #            end coordinate. Default is "end".
  #
  # Returns:
  # A GenomicRanges object for the particular sample that has the breakpoints
  # from the data.frame

  # If no data.frame is provided, stop
  if (is.null(break_df)) {
    stop("No breaks data.frame has been provided. ")
  }

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

  if (sample_id[1] != "all") {
    # Extract samples' info from the data.frame
    break_df <- break_df %>%
      dplyr::filter(!!rlang::sym(samples_col) %in% c(sample_id))

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
    mcols = dplyr::pull(break_df, !!rlang::sym(samples_col))
  )

  return(granges)
}

break_density <- function(breaks_df = NULL,
                          sample_id = NULL,
                          window_size = 1e6,
                          chr_sizes_vector = NULL,
                          samples_col = "samples",
                          chrom_col = "chrom",
                          start_col = "start",
                          end_col = "end") {
  # For given breaks data.frame(s), calculate the breaks density for a tiled
  # windows across the genome. Returns the data as a GenomicRanges object for
  # easy mapping with ggbio. Where the density and counts are stored in
  # @elementMetadata@listData.
  #
  # Args:
  #   breaks_df: a data.frame with chromosomal breaks.
  #   sample_id: a character string that designates which data needs to be
  #              extracted and made intoa GenomicRanges object by matching the
  #              information in the designated sample_col. If "all" is designated,
  #              all samples will be kept. Multiple samples can be designated
  #              as a character vector.
  #   chr_size_vector: a named numeric vector of the sizes (bp) of the chromosomes.
  #                  names of the chromosomes must match the format of the input
  #                  break data. e.g. "chr1" or "1".
  #   window_size: What size windows to calculate break density. Default is 1Mb.
  #   samples_col_sv/cnv: character string that indicates the column name with the
  #                sample ID information. Default is "samples". Will be passed to
  #                `make_granges` function.
  #   chrom_col_sv/cnv: character string that indicates the column name with the
  #              chromosome information. Default is "chrom". Will be passed to
  #              `make_granges`` function.
  #   start_col_sv/cnv: character string that indicates the column name with the
  #              start coordinate. Default is "start". Will be passed to
  #              `make_granges`` function.
  #   end_col_sv/cnv: character string that indicates the column name with the
  #            end coordinate. Default is "end". Will be passed to
  #            `make_granges` function.
  #
  # Check that a sample ID has been specified.
  if (is.null(sample_id)) {
    stop("No sample ID(s) have been specified. Use the `sample_id` argument.")
  }
  # Check that a chromosome sizes have been given
  if (is.null(chr_sizes_vector)) {
    stop("No `chr_sizes_vector` argument has been given. Need that to create bins.")
  }
  # Determine how many samples are in the group
  n_samples <- length(sample_id)

  # Make this into a GenomicRanges object
  break_ranges <- make_granges(
    break_df = breaks_df,
    sample_id = sample_id,
    samples_col = samples_col,
    chrom_col = chrom_col,
    start_col = start_col,
    end_col = end_col
  )
  ######################### Tally breaks by genome bins ########################
  bins <- GenomicRanges::tileGenome(chr_sizes_vector,
    tilewidth = window_size
  )

  # Uncompress GRangesList
  bins <- unlist(bins)

  # Get counts for each genome bin
  bin_counts <- suppressWarnings(GenomicRanges::countOverlaps(
    bins,
    break_ranges
  ))
  ########################### Calculate summary stats ##########################
  # Get a per sample break down if there is more than one sample
  if (n_samples > 1) {

    # Get counts for each genome bin
    bin_indices <- suppressWarnings(GenomicRanges::findOverlaps(
      bins,
      break_ranges
    ))

    # Get list of samples
    bin_samples <- break_ranges@elementMetadata@listData$mcols[bin_indices@to]

    # Make a matrix that has the number of breaks per sample for each bin
    freq_per_bin <- table(bin_indices@from, bin_samples) %>%
      data.frame() %>%
      tidyr::spread(Var1, Freq) %>%
      tibble::column_to_rownames("bin_samples") %>%
      t()

    # Calculate the median breaks per bin
    median_counts <- apply(freq_per_bin, 1, median, na.rm = TRUE)

    # Store the median break counts, some bins data may be dropped off so we need
    # to start with 0s and then fill in the data based on the indices.
    bins@elementMetadata@listData$median_counts <- rep(0, length(bins))
    bins@elementMetadata@listData$median_counts[as.numeric(names(median_counts))] <- median_counts

    # Store average count info
    bins@elementMetadata@listData$avg_counts <- bin_counts / n_samples
  }
  # Store count info
  bins@elementMetadata@listData$total_counts <- bin_counts

  # Calculate and store density
  bins@elementMetadata@listData$density <- bin_counts / window_size

  # Return the GRanges object for mapping purposes
  return(bins)
}
