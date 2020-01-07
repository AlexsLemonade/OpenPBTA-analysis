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

overlap_sv_cnv <- function(cnv_df, 
                           sv_df, 
                           by_sample = TRUE, 
                           max_gap = 0, 
                           samples_col_cnv = "samples", 
                           chrom_col_cnv =  "chrom", 
                           start_col_cnv = "start", 
                           end_col_cnv = "end", 
                           samples_col_sv = "samples", 
                           chrom_col_sv =  "chrom", 
                           start_col_sv = "start", 
                           end_col_sv = "end") {
  # For a given breaks data.frame make a GenomicRanger object from it. Optionally
  # can filter to a single samples' data. 
  #
  # Args:
  #   break_df: for a data.frame with chromosomal coordinates and sample IDs, 
  #             any other columns in this data.frame will also be carried along. 
  #   by_samples: should only samples of the same ID be overlapped? TRUE/FALSE
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
  # Returns:
  # A Genomic Ranges formatted object for the particular sample that has the 
  # break points from the data.frame
  
  # If by_sample is TRUE, obtain common samples list, otherwise make it 
  # null so the whole datasets will be overlapped with no discrimination to 
  # sample. 
  if (!by_sample) {
    common_samples <-  "all"
  } else {
    common_samples <- 
      intersect(dplyr::pull(cnv_df, !!rlang::sym(samples_col_cnv)), 
                dplyr::pull(sv_df, !!rlang::sym(samples_col_sv)))
    
  }
  
  # Obtain a overlap list 
  overlap_granges <- lapply(common_samples, function(sample_id) {
      cnv_range <- make_granges(cnv_df, 
                                sample_id = sample_id, 
                                samples_col = samples_col_cnv,
                                chrom_col =  chrom_col_cnv, 
                                start_col = start_col_cnv, 
                                end_col = end_col_cnv) 
      # Make the SV break data.frame into GRanges objects by each sample
      sv_range <- make_granges(sv_df, 
                               sample_id = sample_id,
                               samples_col = samples_col_sv, 
                               chrom_col =  chrom_col_sv, 
                               start_col = start_col_sv, 
                               end_col = end_col_sv) 
      # Identify overlapping breaks
      overlap <- GenomicAlignments::findOverlaps(cnv_range,
                                                 sv_range, 
                                                 maxgap = max_gap)
      
      # Turn mcols into something useful like what dataset it comes from
      colnames(cnv_range@elementMetadata) <-
        gsub("mcols", "cnv", colnames(cnv_range@elementMetadata))
      
      colnames(sv_range@elementMetadata) <-
        gsub("mcols", "sv", colnames(sv_range@elementMetadata))
      
      # Combine the SV and CNV matches into one data.frame
      break_match <- data.frame(
        cnv_range@elementMetadata[queryHits(overlap), ],
        sv_range@elementMetadata[subjectHits(overlap), ])
      }
  )
  
  # Collapse list into a singular data.frame
  overlap_df <- dplyr::bind_rows(overlap_granges)
  
  return(overlap_df)
}

break_density <- function(sv_breaks = NULL, 
                          cnv_breaks = NULL, 
                          window_size = 1e6,
                          max_gap = 0,
                          chr_sizes_list = NULL,
                          by_sample =
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

map_density_plot <- function(granges = bins, 
                             y_val = "density", 
                             color = "blue", 
                             y_lab = "Breaks per Kb") {
  # Given a GRanges object, plot it y value along its chromosomal mappings using 
  # ggbio. 
  # 
  # Args:
  #   granges: A Granges object to plot
  #   y_val: a character string of the columnname in listData spot of the 
  #          GenomicRanges to plot on the y axis
  #   color: an optional factor level to color code things by
  #   y_lab: a character string to use for the ylabel. Will be passed to 
  #          ggplot2::ylab argument. 
  
  density_plot <- ggbio::autoplot(granges, ggplot2::aes(y = !!rlang::sym(y_val),
                                                        color = !!rlang::sym(color))),
                                  geom = "line", scales = "free_x", space = "free_x") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 3, angle = 45, hjust = 1)) +
    colorblindr::scale_color_OkabeIto(name = color) +
    ggplot2::ylab(y_lab)

  # Print out plot
  density_plot@ggplot
}
