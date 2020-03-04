# Functions for calling CN statuses of genome bins
#
# C. Savonen for ALSF - CCDL
#
# 2020

bp_per_bin <- function(bin_ranges, status_ranges) {
  # Given a binned genome ranges object and another GenomicRanges object, 
  # Return the number of bp covered per bin.
  #
  # Args:
  #   bin_ranges: A binned GenomicRanges made from tileGenome. 
  #   status_ranges:A GenomicRanges object to calculate what percent coverage of
  #   each bin. 
  #
  # Returns:
  #  a data.frame with bins x number of bp 
  
  # Find the portions of each segment that overlap with each bin.
  bin_overlaps <- GenomicRanges::pintersect(
    IRanges::findOverlapPairs(
      bin_ranges,
      status_ranges
    )
  )
  
  # Which segs are a part of which bins?
  bin_indices <- GenomicRanges::findOverlaps(
    bin_overlaps,
    status_ranges
  )
  
  # Get the sum of the length of all excluded regions for each bin.
  bp_per_bin <- tapply(
    bin_overlaps@ranges@width[bin_indices@to], # Get length of each sequence
    bin_indices@from, # Index of which bin it overlaps
    sum
  ) # Add up per bin
  
  # Format as data.frame with rows = bins
  per_bin_df <- data.frame(
    bin = names(bp_per_bin),
    bp_per_bin
  )
  
  # Store dummy counts if there are no ranges that co
  if (nrow(per_bin_df) == 0) {
    per_bin_df <- data.frame(
      bin = factor(1:length(bin_ranges)),
      bp_per_bin = 0
    )
  }
  return(per_bin_df)
}

call_bin_status <- function(sample_id,
                            seg_ranges,
                            bin_ranges,
                            perc_delta_threshold) {
  # Given a sample_id, CN segment ranges, and binned genome ranges object, 
  # make a call for each bin on what CN copy status has the most coverage in the bin. 
  # Uses bp_per_bin function. 
  #
  # Args:
  #   sample_id: A string that corresponds to the biospecimen s
  #   seg_ranges: A GenomicRanges object that contains a `status` and a `biospecimen` slot.  
  #               The `biospecimen slot will be used to split out the `sample_id`'s corresponding ranges.   
  #               The `status` slot should have the that has the gain/loss/neutral 
  #               status call for each segment. Data will be split up by this variable. 
  #   bin_ranges: A binned GenomicRanges made from tileGenome that has been uncompressed with `unlist`. 
  #   perc_delta_threshold: What covereage percent difference do you need to make the call?
  #
  # Returns:
  #  a small data.frame that contains the status call of the sample for each bin. 
  #
  # Extract the ranges for this sample
  sample_seg_ranges <- seg_ranges[which(seg_ranges$biospecimen == sample_id)]
  
  # Split ranges into their respective statuses
  gain_ranges <- sample_seg_ranges[sample_seg_ranges$status == "gain"]
  loss_ranges <- sample_seg_ranges[sample_seg_ranges$status == "loss"]
  neutral_ranges <- sample_seg_ranges[sample_seg_ranges$status == "neutral"]
  
  # Calculate length of each type of status per bin
  gain_per_bin <- bp_per_bin(bin_ranges, gain_ranges)
  loss_per_bin <- bp_per_bin(bin_ranges, loss_ranges)
  neutral_per_bin <- bp_per_bin(bin_ranges, neutral_ranges)
  
  # Format this data into one data.frame where each row is a bin
  bin_bp_status <- data.frame(
    bin = factor(1:length(bin_ranges)),
    # Keep bin width
    bin_width = bin_ranges@ranges@width
  ) %>%
    # Join gains coverage data
    dplyr::left_join(gain_per_bin,
                     by = "bin"
    ) %>%
    # Join loss coverage data
    dplyr::left_join(loss_per_bin,
                     by = "bin",
                     suffix = c(".gain", ".loss")
    ) %>%
    # Join neutral coverage data
    dplyr::left_join(neutral_per_bin,
                     by = "bin"
    ) %>%
    # If there is an NA, at this point we can assume it means 0
    dplyr::mutate_at(
      dplyr::vars(
        dplyr::starts_with("bp_per_bin")
      ),
      ~ tidyr::replace_na(., 0)
    ) %>%
    # Reformat neutral so it is like the others
    dplyr::rename(bp_per_bin.neutral = bp_per_bin) %>%
    # Calculate the bins percentage of each status
    dplyr::mutate(
      perc_gain = bp_per_bin.gain / bin_width,
      perc_loss = bp_per_bin.loss / bin_width,
      perc_neutral = bp_per_bin.neutral / bin_width
    ) %>%
    # Use these percentages for declaring final call per bin based on
    # the perc_delta_threshold
    dplyr::mutate(
      status = dplyr::case_when(
        (perc_gain - perc_loss) > perc_delta_threshold &
          (perc_gain - perc_neutral) > perc_delta_threshold ~ "gain",
        (perc_loss - perc_gain) > perc_delta_threshold &
          (perc_loss - perc_neutral) > perc_delta_threshold ~ "loss",
        TRUE ~ "neutral"
      )
    )
  
  # Format this data as a status
  status_df <- bin_bp_status %>%
    dplyr::select(bin, status) %>%
    tidyr::spread(bin, status) %>%
    dplyr::select(order(as.numeric(colnames(.))))
  
  return(status_df)
}

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
