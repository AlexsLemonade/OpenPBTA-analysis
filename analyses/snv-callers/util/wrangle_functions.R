# Set up and calculate functions for handling MAF data
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################

set_up_maf <- function(maf_df, metadata_df = NULL) {
  # Creates these new variables from a MAF formatted data.frame provided: VAF,
  # mutation_id, base_change, change. Optionally can tack on metadata columns
  # which will be matched using a `Tumor_Sample_Barcode` field. Lastly, any
  # columns that contain all `NA` values will be removed.
  #
  # Args:
  #   maf_df: a maf formatted data.frame
  #   metadata: a data.frame with metadata that you would like to merge with the
  #             maf_df and it's newly calculated variables. (Optional)
  #
  # Returns:
  #   a data.frame with all the original information from the MAF data.frame
  #   maf object but with these new variables: VAF, mutation_id, base_change,
  #   change, coding. If metadata_df is specified, then it will also have
  #   those columns added.

  # Extract the data part of the maf object, put it through a dplyr pipe.
  maf_df <- maf_df %>%
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
    ) %>%
    # Get rid of any variables that have completely NAs.
    dplyr::select(-which(apply(is.na(.), 2, all)))

  # If metadata_df was specified:
  if (!is.null(metadata_df)) {
    # Tack on the metadata so we have this info
    maf_df <- maf_df %>%
      dplyr::left_join(metadata, by = "Tumor_Sample_Barcode") %>%
      # Get rid of any variables that have completely NAs.
      dplyr::select(-which(apply(is.na(.), 2, all)))
  }
}

maf_to_granges <- function(maf_df) {
  # Turn MAF data.frame into a GRanges object. All of the original data.frame will
  # be stored in the `mcols` slot of the GRanges object. This original data
  # frame in the `mcols` slots can be extracted later using
  # @elementMetadata@listData.
  #
  # Args:
  #   maf_df: A MAF formatted data.frame with `Chromosome`, `Start_Position`,
  #           `End_Position`, and maybe `Strand`.
  #   strand: specify if the column `Strand` exists.
  #
  # Returns:
  # A Genomic Ranges formatted object.
  #
  # Create the GRanges object with the strand
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

wxs_bed_filter <- function(maf_df, wxs_bed_file = NULL, bp_window = 0) {
  # Given a MAF formatted data.frame and a BED regions data.frame; filter out
  # any variants of the MAF df that are not within the BED regions.
  #
  # Args:
  #   maf_df: maf data that has been turned into a data.frame. Can be a maf object
  #           that is subsetted using `@data`.
  #   wxs_bed_file: a file path to TSV file that has BED formatted columns with
  #                 chromosome, start, end positions in that order.
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #
  # Returns:
  # The same MAF formatted data.frame with the WXS mutations that lie outside
  # the supplied WXS BED regions filtered out.

  # Read in the BED regions file and make sure the column names are the same
  # as the MAF column format names
  wxs_bed_ranges <- readr::read_tsv(wxs_bed_file, col_names = FALSE) %>%
    dplyr::rename(Chromosome = X1, Start_Position = X2, End_Position = X3)

  # Turn the WXS bed regions into a GRanges object
  wxs_bed_granges <- maf_to_granges(wxs_bed_ranges)

  # Obtain a MAF data.frame of only the WXS samples since this filter will only
  # be applied to those samples.
  maf_wxs <- maf_df %>%
    dplyr::filter(experimental_strategy == "WXS")

  # Error catcher in case there are no `WXS` samples
  if (nrow(maf_wxs) == 0) {
    stop("No WXS samples found underneath column 'experimental_strategy'
            double check filtering steps and data.")
  }

  # Turn the MAF WXS sample mutations into a GRanges object
  wxs_maf_granges <- maf_to_granges(maf_wxs)

  # Find the overlap of the BED regions and the mutations This outputs a
  # special GenomicRanges object that contains indices of each of these
  # ranges that overlap
  overlap <- GenomicRanges::findOverlaps(
    wxs_maf_granges,
    wxs_bed_granges,
    maxgap = bp_window
  )

  # Calculate of ratio of variants in this BED using the @from slot which
  # indicates the indices of the ranges in `wxs_maf_ranges` that have overlaps
  # with `wxs_bed_ranges`
  ratio <- length(overlap@from) / nrow(maf_wxs)

  # What fraction of mutations are in these bed regions?
  cat(
    "Ratio of variants in this BED:", ratio, "\n",
    "Ratio of variants being filtered out:", 1 - ratio
  )

  # Only keep those in the BED regions that overlap the `wxs_bed_granges`
  maf_wxs <- maf_wxs[unique(overlap@from), ]

  # Tack these back on to the WGS samples which remain unfiltered
  filt_maf_df <- maf_df %>%
    dplyr::filter(experimental_strategy == "WGS") %>%
    dplyr::bind_rows(maf_wxs)

  # Return this matrix with the WXS mutations filtered but WGS the same
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
  #   maf_df: maf data.frame that has been turned into a data.frame, has had
  #           the experimental_stategy column added from the metadata (can be
  #           done with the `set_up_maf` function) and has WXS mutations filtered
  #           using `wxs_bed_filter` function (If the situation calls for it).
  #   wgs_size: genome size in bp to be used for WGS samples
  #   wxs_size: genome size in bp to be used for WGS samples
  #
  # Returns:
  # A sample-wise data.frame with Tumor Mutation Burden statistics calculated
  # using the given WGS and WXS sizes.
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

annotr_maf <- function(maf_df, annotation_file = NULL, bp_window = 0) {
  # Annotate the genomic regions mutations are in from a MAF formatted data.
  # frame with AnnotatR object. Because any given region can have multiple
  # region labels, a single mutation in a can have a single label, multiple
  # labels or no labels. Generally this data.frame ends up being very large
  # because many mutations have several labels.
  #
  # Args:
  #   maf_df: maf data that has been turned into a data.frame.
  #   annotation_file: a file path to a .RDS object containing the desired annotation
  #                    from AnnotatR build function.
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #
  # Returns:
  # A large data.frame that contains every mutation and every region type it
  # overlaps. Genomic region types are noted in the `type` column.
  #
  # Read in the genomic regions annotation object
  annotation_ranges <- readr::read_rds(annotation_file)

  # Use custom function to turn our MAF data into a GRanges
  maf_ranges <- maf_to_granges(maf_df)

  # Intersect the regions we read in with the annotations.
  # This outputs a special GenomicRanges object that contains indices of
  # each of these ranges that overlap
  annot_matches <- GenomicRanges::findOverlaps(
    annotation_ranges,
    maf_ranges,
    maxgap = bp_window
  )

  # Extract the annotation from the regions that overlap the mutations
  annot <- data.frame(
    annotation_ranges[annot_matches@from]@elementMetadata@listData, # Extract matching annotation ranges info
    maf_ranges[annot_matches@to]@elementMetadata@listData, # Extract matching MAF ranges data.frame info
    stringsAsFactors = FALSE
  ) %>%
    # We will rename all columns that start with `mcols` and drop the `mcols` string
    # Because `mcols.` has 6 characters, we start at 7, and want to make to get all the
    # characters in the original name, so we just say a number that is much bigger than
    # the longest column name (like 100)
    dplyr::rename_at(dplyr::vars(dplyr::starts_with("mcols.")), substr, 7, 100) %>%
    dplyr::mutate("type" = as.factor(gsub("^hg38_genes_", "", type)))

  return(annot)
}

find_cosmic_overlap <- function(maf_df, cosmic_clean_file, bp_window = 0) {
  # Find how much overlap a MAF formatted data.frame has with COSMIC mutations
  # Return columns in the MAF data.frame that say whether or not it is also
  # in the COSMIC mutations set.
  #
  # Args:
  #   maf_df: MAF formatted data that has been turned into a data.frame and has
  #           been run through `set_up_variables`
  #   cosmic_clean_file: a file path to a TSV of COSMIC mutations that has been
  #                     been cleaned up to have the genomic coordinates separated
  #                     into Chr, Start, and End columns.
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #
  # Returns:
  # The original MAF data.frame is returned with two added columns, called
  # overlap_cosmic and `same_cosmic`. These columns are logical type and
  # indicate whether the mutation was overlapping a COSMIC mutation
  # (`overlap_cosmic`) and whether it was overlapping a COSMIC mutation and also
  # contained the same overall base change (`same_cosmic`).
  
  # Read in the data
  cosmic_df <- readr::read_tsv(cosmic_clean_file)
  
  # Identify which rows have NAs
  cosmic_na <- apply(is.na(cosmic_df[, 1:4]), 1, sum)
  
  # Keep only the mutations with adequate information about their locations
  cosmic_df <- cosmic_df %>% 
    dplyr::filter(cosmic_na < 1)
  
  # Read in the cosmic file and turn it into a GRanges object
  cosmic_granges <- maf_to_granges(cosmic_df)

  # Turn the MAF data.frame into a GRanges object
  maf_granges <- maf_to_granges(maf_df)

  # Find the overlap of the MAF and COSMIC mutations.
  # This outputs a special GenomicRanges object that contains indices of
  # each of these ranges that overlap
  overlap <- GenomicRanges::findOverlaps(
    maf_granges,
    cosmic_granges,
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
