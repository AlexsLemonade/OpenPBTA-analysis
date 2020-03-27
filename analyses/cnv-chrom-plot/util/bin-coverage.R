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
  #   status_ranges: A GenomicRanges object to calculate what percent coverage of
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
    bin_ranges,
    bin_overlaps
  )
  
  # Get the sum of the length of all seg portions for each bin.
  bp_per_bin <- tapply(
    bin_overlaps@ranges@width, # Get length of each sequence within the bin
    bin_indices@from, # Index of which bin it overlaps
    sum
  ) # Add up length per bin
  
  # Format as data.frame with rows = bins
  per_bin_df <- data.frame(
    bin = as.numeric(names(bp_per_bin)),
    bp_per_bin = as.numeric(bp_per_bin)
  )
  
  # Store dummy counts if there are no ranges that are in the bins
  if (nrow(per_bin_df) == 0) {
    per_bin_df <- data.frame(
      bin = as.numeric(1:length(bin_ranges)),
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
    bin = as.numeric(1:length(bin_ranges)),
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
        (perc_gain - perc_loss) > perc_delta_threshold ~ "gain",
        (perc_loss - perc_gain) > perc_delta_threshold ~ "loss",
        TRUE ~ "neutral"
      )
    )
  
  # Format this data as a status
  status_df <- bin_bp_status %>%
    # Only keep the bin and status columns
    dplyr::select(bin, status) %>%
    # Spread this data so we can make it a sample x bin matrix later
    tidyr::spread(bin, status) %>%
    # Sort the bins data/columns to be in numeric order
    dplyr::select(order(as.numeric(colnames(.))))
  
  return(status_df)
}
