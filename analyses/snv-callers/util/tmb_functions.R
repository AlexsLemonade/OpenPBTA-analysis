# Set up and calculate functions for handling MAF data
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################

maf_to_granges <- function(maf_df) {
  # Turn MAF data.frame into a GRanges object. All of the original data.frame will
  # be stored in the `mcols` slot of the GRanges object. This original data
  # frame in the `mcols` slots can be extracted later using
  # @elementMetadata@listData.
  #
  # Args:
  #   maf_df: A MAF formatted data.frame with `Chromosome`, `Start_Position`,
  #           `End_Position`
  #
  # Returns:
  # A Genomic Ranges formatted object.
  #
  # Create the GRanges object
  GenomicRanges::GRanges(
    seqnames = maf_df$Chromosome,
    ranges = IRanges::IRanges(
      start = maf_df$Start_Position,
      end = maf_df$End_Position
    ),
    mcols = maf_df
  )
}

snv_ranges_filter <- function(maf_df, keep_ranges = NULL, bp_window = 0) {
  # Given a MAF formatted data.frame and a BED regions data.frame; filter out
  # any variants of the MAF df that are not within the BED regions.
  #
  # Args:
  #   maf_df: maf data that has been turned into a data.frame. Can be a maf object
  #           that is subsetted using `@data`.
  #   keep_ranges: BED ranges data.frame with columns: chromosome, start, end
  #             positions in that order or a Genomic Ranges object. If data.frame,
  #             is given, it will be converted to GRanges object
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #
  # Returns:
  # The same MAF formatted data.frame with the mutations that lie outside
  # the supplied BED regions filtered out.

  # Turn the MAF sample mutations into a GRanges object
  maf_granges <- maf_to_granges(maf_df)
  
  # If ranges is given as a data.frame, convert
  if (is.data.frame(keep_ranges)) {
    # Turn the bed regions df into a GRanges object
    keep_ranges <- GenomicRanges::GRanges(
      seqnames = keep_ranges$X1,
      ranges = IRanges::IRanges(
        start = keep_ranges$X2,
        end = keep_ranges$X3
      )
    )
  }
  # Find the overlap of the BED regions and the mutations This outputs a
  # special GenomicRanges object that contains indices of each of these
  # ranges that overlap
  overlap <- GenomicRanges::findOverlaps(
    maf_granges,
    keep_ranges,
    maxgap = bp_window
  )
  
  # Calculate of ratio of variants in this BED using the @from slot which
  # indicates the indices of the ranges in `wxs_maf_ranges` that have overlaps
  # with `wxs_bed_ranges`
  ratio <- length(overlap@from) / nrow(maf_df)
  
  # What fraction of mutations are in these bed regions?
  cat(
    "Ratio of variants in this BED:", ratio, "\n",
    "Ratio of variants being filtered out:", 1 - ratio, "\n"
  )
  
  # Only keep those in the BED regions that overlap the `wxs_bed_granges`
  filt_maf_df <- maf_df[unique(overlap@from), ]
  
  # Return this matrix with the WXS mutations filtered but WGS the same
  return(filt_maf_df)
}

calculate_tmb <- function(maf_df, bed_wgs, bed_wxs) {
  # Calculate Tumor Mutational Burden for each sample given their WGS or WXS
  # sizes in bp. Based on the BED files provided, tthis function filters out SNV 
  # outside those ranges and uses the total size in bp of those BED ranges for th
  # TMB denominator. 
  #
  # TMB = # variants / size of the genome or exome surveyed
  #
  # Args:
  #   maf_df: maf data.frame that has been turned into a data.frame, has had
  #           the experimental_stategy column added from the metadata (can be
  #           done with the `set_up_maf` function) and has WXS mutations filtered
  #           using `wxs_bed_filter` function (If the situation calls for it).
  #   bed_wgs: BED file path with genome ranges to be used for filtering and 
  #            genome size denominator for WGS samples
  #   bed_wxs: BED file path genome ranges to be used for filtering and genome 
  #            size denominator for WXS samples
  #
  # Returns:
  # A sample-wise data.frame with Tumor Mutation Burden statistics calculated
  # using the given WGS and WXS sizes.
  
  # Read in the BED files we need
  bed_wgs <- readr::read_tsv(bed_wgs, col_names = FALSE)
  bed_wxs <- readr::read_tsv(bed_wxs, col_names = FALSE)
  
  # Sum up genome sizes
  wgs_size <- sum(bed_wgs[, 3] - bed_wgs[, 2])
  wxs_size <- sum(bed_wxs[, 3] - bed_wxs[, 2])
  
  # Don't want integers per se
  wgs_size <- as.numeric(wgs_size)
  wxs_size <- as.numeric(wxs_size)
  
  # For WXS samples, filter out mutations that are outside of these coding regions.
  filt_wxs_maf_df <- snv_ranges_filter(dplyr::filter(
    maf_df,
    experimental_strategy == "WXS"
  ),
  keep_ranges = bed_wxs
  )
  
  # For WXS samples, filter out mutations that are outside of these coding regions.
  filt_wgs_maf_df <- snv_ranges_filter(dplyr::filter(
    maf_df,
    experimental_strategy == "WGS"
  ),
  keep_ranges = bed_wgs
  )
  
  # Bind the filtered WXS sample rows back to the WGS samples
  filt_maf_df <- dplyr::bind_rows(filt_wxs_maf_df, filt_wgs_maf_df)
  
  # Make a genome size variable
  tmb <- maf_df %>%
    dplyr::mutate(genome_size = dplyr::recode(experimental_strategy,
                                              "WGS" = wgs_size,
                                              "WXS" = wxs_size
    )) %>%
    dplyr::group_by(
      Tumor_Sample_Barcode, experimental_strategy, genome_size,
      short_histology
    ) %>%
    # Count number of mutations for that sample
    dplyr::summarize(mutation_count = dplyr::n()) %>%
    # Calculate TMB
    dplyr::mutate(tmb = mutation_count / (genome_size / 1000000))
  
  return(tmb)
}
