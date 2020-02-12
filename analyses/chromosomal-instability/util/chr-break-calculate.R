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
                          end_col = "end",
                          unsurveyed_bed = NULL,
                          perc_cutoff = .75,
                          return_vector = FALSE) {
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
  #   unsurveyed_bed: an optional BED format data frame with columns `chrom`, `start`
  #               and `end` columns that specify regions that should be NA
  #               because they were not surveyed.
  #   perc_cutoff: Only relevant of `unsurveyed_bed` is being used. The percent
  #                overlap above which a bin's region should be considered NA.
  #   return_vector: Default is to return a GenomicRanges, but if this option is
  #                  used, only a named vector is returned with the data for
  #                  each bin
  # Returns:
  # Breakpoint density and total counts per bin are stored either as a GenomicRanges
  # object(default) @elementMetadata@listData or a named vector.
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
  bins@elementMetadata@listData$density <- bin_counts / bins@ranges@width

  # If specified, make unsurveyed into a GenomicRanges
  if (!is.null(unsurveyed_bed)) {
    # Make GRanges for uncallable ranges
    unsurveyed_ranges <- GenomicRanges::GRanges(
      seqnames = unsurveyed_bed$chrom,
      ranges = IRanges::IRanges(
        start = unsurveyed_bed$start,
        end = unsurveyed_bed$end
      )
    )

    # Find the NA region(s) for each bin without combining close regions
    na_regions <- GenomicRanges::pintersect(IRanges::findOverlapPairs(bins, unsurveyed_ranges))

    # Find overlap between na_regions and bin
    na_overlaps <- GenomicRanges::findOverlaps(bins, na_regions)

    # Get the sum of the length of all excluded regions for each bin.
    excluded_length_per_bin <- tapply(
      na_regions@ranges@width, # Get length of each sequence
      na_overlaps@from, # Index of which bin it overlaps
      sum
    ) # Add up per bin

    # Get the total bin length for each bin that has excluded regions
    bin_length <- bins[unique(na_overlaps@from)]@ranges@width

    # Calculate the percent overlap, adding up the na regions within a bin
    pct_overlap <- excluded_length_per_bin / bin_length

    # Store bins as names for sanity checking
    # Note that unique and tapply put the bins in the same order
    names(pct_overlap) <- unique(na_overlaps@from)

    # Get the bin indices that correspond to less than the cutoff
    bin_indices <- names(pct_overlap)[which(pct_overlap > perc_cutoff)]

    # Make sure these are numeric again so we can reference them
    bin_indices <- as.numeric(bin_indices)

    # Call these data points NA instead
    bins$total_counts[bin_indices] <- NA
    bins$density[bin_indices] <- NA
  }

  # If return_vector option is used, extract only the data
  if (return_vector) {
    # Just extract vector of data
    total_counts <- bins@elementMetadata@listData$total_counts

    # Name with chromosome labels
    names(total_counts) <- S4Vectors::decode(bins@seqnames)

    return(total_counts)
  } else {
    # Return the GRanges object for mapping purposes
    return(bins)
  }
}
