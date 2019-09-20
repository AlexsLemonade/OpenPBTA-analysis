# Set up and calculate functions for handling MAF data
#
# C. Savonen for ALSF - CCDL
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################

calculate_vaf <- function(maf_df) {
  # Creates these new variables from a MAF formatted data.frame provided: VAF,
  # mutation_id, base_change, change.
  #
  # Args:
  #   maf_df: a maf formatted data.frame
  #
  # Returns:
  #   a data.frame with all the original information in the `@data` part of the
  #   maf object but with these new variables: VAF, mutation_id, base_change,
  #   change, coding.

  # Extract the data part of the maf object, put it through a dplyr pipe.
  maf_df %>%
    dplyr::mutate(
      # Calculate the variant allele frequency
      vaf = as.numeric(t_alt_count) / (as.numeric(t_ref_count) +
        as.numeric(t_alt_count)),
      # Create a base_change variable
      base_change = paste0(Reference_Allele, ">", Allele),

      # Create a numeric portion of the PolyPhen score
      PolyPhen_numeric = as.numeric(stringr::word(PolyPhen, 2, sep = "\\(|\\)")),

      # Create a categorical portion of the PolyPhen score
      PolyPhen_category = stringr::word(PolyPhen, 1, sep = "\\("),

      Variant_Classification = as.factor(Variant_Classification)
    ) %>%
    dplyr::mutate(
      # From the base_change variable, summarize insertions, deletions, and
      # changes that are more than one base into their own groups.
      change = dplyr::case_when(
        grepl("^-", base_change) ~ "ins",
        grepl("-$", base_change) ~ "del",
        nchar(base_change) > 3 ~ "long_change",
        TRUE ~ base_change
      )
    ) %>%
    dplyr::mutate(
      # Create the mutation id based on the change variable as well as the
      # gene symbol, start position, and sample ID.
      mutation_id = paste0(
        Hugo_Symbol, "_",
        change, "_",
        Start_Position, "_",
        Tumor_Sample_Barcode
      ),
    )
}

maf_to_granges <- function(maf_df) {
  # Turn MAF data.frame into a GRanges object. All of the original data.frame will
  # be stored in the `mcols` slot of the GRanges object.
  #
  # Args:
  #   maf_df: A MAF formatted data.frame with `Chromosome`, `Start_Position`,
  #           `End_Position`, and `Strand`.
  #
  # Create the GRanges object
  GenomicRanges::GRanges(
    seqnames = maf_df$Chromosome,
    ranges = IRanges::IRanges(
      start = maf_df$Start_Position,
      end = maf_df$End_Position
    ),
    strand = maf_df$Strand,
    mcols = maf_df
  )
}

wxs_bed_filter <- function(maf_df, wxs_bed = NULL, bp_window = 0) {
  # Given a MAF formatted data.frame and a BED regions data.frame; filter out
  # any variants of the MAF df that are not within the BED regions.
  #
  # Args:
  #   maf_df: maf data that has been turned into a data.frame. Can be a maf object
  #           that is subsetted using `@data`.
  #   wxs_bed: a data.frame that has the windows to be used for TMB calculation.
  #            BED formatted columns with chromosome, start, end positions in
  #            that order.
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.

  # Turn the WXS bed regions into a GRanges object
  wxs_bed_granges <- GenomicRanges::GRanges(
    seqnames = wxs_bed$X1,
    ranges = IRanges::IRanges(
      start = wxs_bed$X2,
      end = wxs_bed$X3
    )
  )

  # Obtain a MAF data.frame of only the WXS samples
  maf_wxs <- maf_df %>%
    dplyr::filter(experimental_strategy == "WXS")

  # Error catcher
  if (nrow(maf_wxs) == 0) {
    warning("No WXS samples found underneath column 'experimental_strategy'
            double check filtering steps and data.")
  }

  # Find the overlap of the BED regions and the mutations.
  wxs_maf_granges <- GenomicRanges::GRanges(
    seqnames = maf_wxs$Chromosome,
    ranges = IRanges::IRanges(
      start = maf_wxs$Start_Position,
      end = maf_wxs$End_Position
    )
  )

  # Find the overlap of the BED regions and the mutations.
  overlap <- GenomicRanges::findOverlaps(wxs_maf_granges, wxs_bed_granges,
    maxgap = bp_window
  )

  # Calculate of ratio of variants in this BED
  ratio <- length(overlap@from) / nrow(maf_wxs)

  # What fraction of mutations are in these bed regions?
  cat(
    "Ratio of variants in this BED:", ratio,
    "Ratio of variants being filtered out:", 1 - ratio
  )

  # Only keep those in the BED regions
  maf_wxs <- maf_wxs[unique(overlap@from), ]

  # Tack these back on to the WGS samples
  filt_maf_df <- maf_df %>%
    dplyr::filter(experimental_strategy == "WGS") %>%
    dplyr::bind_rows(maf_wxs)

  # Return this filtered matrix
  return(filt_maf_df)
}

calculate_tmb <- function(maf_df, wgs_size, wxs_size) {
  # Calculate Tumor Mutational Burden for each sample given their WGS or WXS
  # sizes in bp. Ideally have filtered the WXS mutations previously using the
  # `wxs_bed_filter` function.
  #
  # TMB = # variants / size of the genome or exome surveyed
  #
  # Args:
  #   maf_df: maf data.frame that has been turned into a data.frame and has
  #           WXS mutations filtered using `wxs_bed_filter` function.
  #   wgs_size: genome size in bp to be used for WGS samples
  #   wxs_size: genome size in bp to be used for WGS samples
  #
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

annotr_maf <- function(maf_df, annotation_obj = annotations, bp_window = 0) {
  # Annotate the genomic regions mutations are in from a MAF formatted data.
  # frame with AnnotatR object.
  #
  # Args:
  #   maf_df: maf data that has been turned into a data.frame.
  #   annotation_obj: an annotation object with the desired annotations.
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #
  # Read in the genomic regions annotation object
  annotations <- readr::read_rds(annot_rds)

  # Use custom function to turn our MAF data into a GRanges
  maf_ranges <- maf_to_granges(maf_df)

  # Intersect the regions we read in with the annotations
  annot_matches <- GenomicRanges::findOverlaps(annotation_obj, maf_ranges,
    maxgap = bp_window
  )

  # Extract the annotation from the regions that overlap the mut
  annot <- data.frame(
    annotation_obj[annot_matches@from]@elementMetadata@listData,
    maf_ranges[annot_matches@to]@elementMetadata@listData,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::rename_at(dplyr::vars(dplyr::starts_with("mcols.")), substr, 7, 10000) %>%
    dplyr::mutate("type" = as.factor(gsub("^hg38_genes_", "", type)))

  return(annot)
}

find_cosmic_overlap <- function(maf_df, bp_window = 0) {
  # Find how much overlap a MAF formatted data.frame has with COSMIC mutations
  # Return a column in the MAF data.frame that says whether or not it is also
  # in the COSMIC mutations set.
  #
  # Args:
  #   maf_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `set_up_variables`
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #
  # Turn both these datasets into GRanges
  maf_granges <- maf_to_granges(maf_df)
  cosmic_granges <- maf_to_granges(suppressMessages(readr::read_tsv(cosmic_clean_file)))

  # Find the overlap of the MAF and COSMIC mutations.
  overlap <- GenomicRanges::findOverlaps(maf_granges, cosmic_granges,
    maxgap = bp_window
  )

  # Check if the base change is the same
  same_change <- maf_granges@elementMetadata@listData$mcols.change[overlap@from] ==
    cosmic_granges@elementMetadata@listData$mcols.base_change[overlap@to]

  # Make a list of mutation ids that overlap COSMIC mutations
  overlap_w_cosmic <- maf_granges@elementMetadata@listData$mcols.mutation_id[overlap@from]

  #  Make a list of mutation ids that also have the same base_change
  same_as_cosmic <- overlap_w_cosmic[which(same_change)]

  # Calculate of ratio of mutations that overlap with COSMIC mutations
  ratio <- length(unique(overlap_w_cosmic)) / nrow(maf_df)

  # What fraction of mutations are in these bed regions?
  cat(
    " Ratio of variants overlapping with COSMIC:", ratio, "\n",
    "Number of mutations with same base_change:", sum(same_change)
  )

  # Make this a new column
  maf_df %>%
    dplyr::mutate(
      overlap_cosmic = (mutation_id %in% overlap_w_cosmic),
      same_cosmic = (mutation_id %in% same_as_cosmic)
    )
}
