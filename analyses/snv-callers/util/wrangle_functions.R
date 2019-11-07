# Set up and calculate functions for handling MAF data
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Import special functions
source(file.path(root_dir, "analyses", "snv-callers", "util", "sql_functions.R"))

set_up_maf <- function(maf_path = opt$maf, 
                       sql_file = opt$sql_file, 
                       label = opt$label, 
                       overwrite_it = opt$overwrite,
                       vaf_cutoff = 0) {
  # Calculates and filters by VAF
  #
  # Args:
  #   maf_path: a maf formatted data.frame.
  #   sql_file: path to sql file that the metadata is saved to. 
  #   label: what the dataset should be called in the SQL file.
  #   overwrite_it: should it be overwritten in the SQL file? 
  #   vaf_cutoff: what is the minimum VAF that should be kept?
  #
  # Returns:
  #   a data.frame with all the original information from the MAF data.frame
  #   maf object but with these new variables: VAF, mutation_id, base_change,
  #   change, coding. If metadata_df is specified, then it will also have
  #   those columns added.
  #   
  # Save the original 
  save_to_sql(data.table::fread(file.path(opt$maf), skip = 1), 
              sql_file = sql_file, 
              tbl_name = paste0(opt$label, "_vaf"))
  
  # Establish connection
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_file)
  
  # Let's make the SQL table its own object
  maf_df_db <- dplyr::tbl(con, paste0(label, "_vaf"))
  
  # We will do the same for the metadata
  # This should have been set up in 00-set_up.R
  metadata_db <- dplyr::tbl(con, "metadata")
  
  # Extract the data part of the maf object, put it through a dplyr pipe.
  maf_df_db <- maf_df_db %>%
    # Make these numeric
    dplyr::mutate(t_alt_count = as.numeric(t_alt_count), 
                  t_ref_count = as.numeric(t_ref_count), 
                  t_alt_count = as.numeric(t_alt_count)
                  ) %>%
    # Calculate the variant allele frequency
    dplyr::mutate(
      vaf = t_alt_count / 
        (t_ref_count + t_alt_count)) %>% 
    
    # Make the base_change variable
    dplyr::mutate(
      base_change = paste0(Reference_Allele, ">", Allele)) %>% 
    
    # Filter by VAF
    dplyr::filter(vaf > vaf_cutoff, 
                 !is.na(vaf)) %>%
      dplyr::left_join(metadata_db, by = "Tumor_Sample_Barcode") 

  # Collect the resulting data as a data.frame and return it back
  maf_df <- maf_df_db %>% 
    dplyr::collect() %>% 
    
    # This last step can't be done in SQL, so needs to be done in the R env
    dplyr::mutate(
      # From the base_change variable, summarize insertions, deletions, and
      # changes that are more than one base into their own groups.
      change = dplyr::case_when(
        grepl("^-", base_change) ~ "ins",
        grepl("-$", base_change) ~ "del",
        nchar(base_change) > 3 ~ "long_change",
        TRUE ~ base_change
      )) 
  
  # Save the original 
  save_to_sql(maf_df, 
              sql_file = sql_file, 
              tbl_name = paste0(opt$label, "_vaf"))
  
  # Return the data.frame
  message(paste(nrow(maf_df), "mutations left after filter and merge"))
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
  #
  # Returns:
  # A Genomic Ranges formatted object.

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

wxs_bed_filter <- function(sql_file = opt$sql_file, 
                           label = opt$label, 
                           overwrite_it = opt$overwrite,
                           wxs_bed_file = opt$bed_wxs, 
                           bp_window = 0) {
  # Given a MAF formatted data.frame and a BED regions data.frame; filter out
  # any variants of the MAF df that are not within the BED regions.
  #
  # Args:
  #   sql_file: path to sql file that the metadata is saved to. 
  #   label: what the dataset should be called in the SQL file.
  #   wxs_bed_file: a file path to TSV file that has BED formatted columns with
  #                 chromosome, start, end positions in that order.
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #   overwrite_it: should it be overwritten in the SQL file? 
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

  # Extract the df from SQL file
  maf_df <- sql_to_df(sql_file = sql_file, 
                      tbl_name = paste0(opt$label, "_vaf"))

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
    "Ratio of variants being filtered out:", 1 - ratio, "\n"
  )

  # Only keep those in the BED regions that overlap the `wxs_bed_granges`
  maf_wxs <- maf_wxs[unique(overlap@from), ]

  # Tack these back on to the WGS samples which remain unfiltered
  filt_maf_df <- maf_df %>%
    dplyr::filter(experimental_strategy == "WGS") %>%
    dplyr::bind_rows(maf_wxs)
 
  # Save this to the SQL file
  save_to_sql(filt_maf_df, 
              sql_file = sql_file, 
              tbl_name = paste0(opt$label, "_vaf"))
}

calculate_tmb <- function(sql_file = opt$sql_file, 
                          label = opt$label,
                          wgs_size, 
                          wxs_size, 
                          overwrite_it = opt$overwrite) {
  # Calculate Tumor Mutational Burden for each sample given their WGS or WXS
  # sizes in bp. Ideally have filtered the WXS mutations previously using the
  # `wxs_bed_filter` function.
  #
  # TMB = # variants / size of the genome or exome surveyed
  #
  # Args:
  #   sql_file: path to sql file that the metadata is saved to. 
  #   label: what the dataset should be called in the SQL file.
  #   wgs_size: genome size in bp to be used for WGS samples
  #   wxs_size: genome size in bp to be used for WGS samples
  #   overwrite_it: should it be overwritten in the SQL file? 
  #
  # Returns:
  # A sample-wise data.frame with Tumor Mutation Burden statistics calculated
  # using the given WGS and WXS sizes.
  # Save to SQL Lite file
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_file)
  
  # Let's make the SQL table its own object
  maf_df_db <- dplyr::tbl(con, paste0(label, "_vaf"))
  
  # Make a genome size variable
  tmb_db <- maf_df_db %>%
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
 
  # Save this to the SQL file
  save_to_sql(tmb_db, 
              sql_file = sql_file, 
              tbl_name = paste0(opt$label, "_tmb"))
}

annotr_maf <- function(annotation_file = opt$annot_rds, 
                       sql_file = opt$sql_file, 
                       label = opt$label, 
                       overwrite_it = opt$overwrite, 
                       bp_window = 0) {
  # Annotate the genomic regions mutations are in from a MAF formatted data.
  # frame with AnnotatR object. Because any given region can have multiple
  # region labels, a single mutation in a can have a single label, multiple
  # labels or no labels. Generally this data.frame ends up being very large
  # because many mutations have several labels.
  #
  # Args:
  #   sql_file: path to sql file that the metadata is saved to. 
  #   label: what the dataset should be called in the SQL file.
  #   annotation_file: a file path to a .RDS object containing the desired annotation
  #                    from AnnotatR build function.
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #   overwrite_it: should it be overwritten in the SQL file? 
  #
  # Returns:
  # A large data.frame that contains every mutation and every region type it
  # overlaps. Genomic region types are noted in the `type` column.

  # Read in the genomic regions annotation object
  annotation_ranges <- readr::read_rds(annotation_file)

  # Read in the maf_df
  maf_df <- sql_to_df(sql_file = sql_file, 
                      tbl_name = paste0(opt$label, "_vaf"))
  
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
    dplyr::mutate("type" = as.factor(gsub("^hg38_genes_", "", type))) %>% 
    # SQL doesn't like having two symbol columns even though one is capitalized and the 
    # other isn't. 
    dplyr::select(-SYMBOL)
   
  # Save this to the SQL file
  save_to_sql(annot, 
              sql_file = sql_file, 
              tbl_name = paste0(opt$label, "_annot"))
  
  # Remove the annotation ranges file to conserve memory burden
  rm(annotation_ranges)
  rm(annot)
  
  # Return annotated mutations
  return(annot)
}

find_cosmic_overlap <- function(cosmic_clean_file, 
                                sql_file = opt$sql_file,
                                label = opt$label, 
                                overwrite_it = opt$overwrite,
                                bp_window = 0) {
  # Find how much overlap a MAF formatted data.frame has with COSMIC mutations
  # Return columns in the MAF data.frame that say whether or not it is also
  # in the COSMIC mutations set.
  #
  # Args:
  #   cosmic_clean_file: a file path to a TSV of COSMIC mutations that has been
  #                     been cleaned up to have the genomic coordinates separated
  #                     into Chr, Start, and End columns.
  #   sql_file: path to sql file that the metadata is saved to. 
  #   label: what the dataset should be called in the SQL file.
  #   overwrite_it: should it be overwritten in the SQL file? 
  #   vaf_cutoff: what is the minimum VAF that should be kept?
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

  # Read in the maf_df
  maf_df <- sql_to_df(sql_file = sql_file, 
                      tbl_name = paste0(opt$label, "_vaf"))
  
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
  overlap_w_cosmic <-
    maf_granges@elementMetadata@listData$mcols.mutation_id[overlap@from]

  #  Make a list of mutation ids that also have the same base_change
  same_as_cosmic <- overlap_w_cosmic[which(same_change)]

  # Calculate of ratio of mutations that overlap with COSMIC mutations
  ratio <- length(unique(overlap_w_cosmic)) / nrow(maf_df)

  # What fraction of mutations are in these bed regions?
  cat(
    " Ratio of variants overlapping with COSMIC:", ratio, "\n",
    "Number of mutations with same base_change:", sum(same_change), "\n"
  )

  # Remove these files to reduce memory burden
  rm(cosmic_df)
  rm(cosmic_na)

  # Make this a new column
  maf_df %>%
    dplyr::mutate(
      overlap_cosmic = (mutation_id %in% overlap_w_cosmic),
      same_cosmic = (mutation_id %in% same_as_cosmic)
    )
  
  # Save this to the SQL file
  save_to_sql(maf_df, 
              sql_file = sql_file, 
              tbl_name = paste0(opt$label, "_vaf"))
}
