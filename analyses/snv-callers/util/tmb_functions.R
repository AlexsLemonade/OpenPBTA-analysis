# Functions for calculating tumor mutation burden
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

snv_ranges_filter <- function(maf_df,
                              keep_ranges = NULL,
                              bp_window = 0) {
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
  # indicates the indices of the ranges in `maf_ranges` that have overlaps
  # with `keep_ranges`
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

calculate_tmb <- function(tumor_sample_barcode = NULL,
                          maf_df,
                          bed_ranges) {
  # Calculate Tumor Mutational Burden a given sample in `Tumor_Sample_Barcode`
  # given a target region BED frame. This function uses `snv_ranges_filter` to
  # filter out SNVs outside those target BED ranges and uses the total size in
  # bp of those BED ranges for the TMB denominator.
  #
  # TMB = # variants / size of the genome or exome surveyed / Mb
  #
  # Args:
  #   tumor_sample_barcode: A string corresponding to a sample
  #                         id in the `Tumor_Sample_Barcode` MAF file in maf_df.
  #   maf_df: maf data.frame that has been turned into a data.frame, has had
  #           the experimental_stategy column added from the metadata (can be
  #           done with the `set_up_maf` function) and has WXS mutations filtered
  #           using `wxs_bed_filter` function (If the situation calls for it).
  #   bed_ranges: A GenomicRanges made from the BED file of the target regions to use for TMB
  #           calculations.
  #
  # Returns:
  # A calculated total region size, and TMB for the given Tumor_Sample_Barcode,
  # returned as a single row data.frame. `experimental_strategy`, `short_histology`
  # columns are carried along.

  # Sum up genome sizes
  bed_size <- as.numeric(sum(bed_ranges@ranges@width))

  # Filter to only the sample's mutations
  sample_maf_df <- maf_df %>%
    dplyr::filter(Tumor_Sample_Barcode == tumor_sample_barcode)

  # Filter out mutations that are outside of these coding regions.
  filt_maf_df <- snv_ranges_filter(sample_maf_df, keep_ranges = bed_ranges)

  # Make a genome size variable
  tmb <- filt_maf_df %>%
    dplyr::group_by(
      #TODO: Make this column passing stuff more flexible with some tidyeval maybe
      Tumor_Sample_Barcode = tumor_sample_barcode,
      experimental_strategy,
      short_histology
    ) %>%
    # Count number of mutations for that sample
    dplyr::summarize(
      mutation_count = dplyr::n(),
      region_size = bed_size,
      tmb = mutation_count / (region_size / 1000000)
      )


  return(tmb)
}
