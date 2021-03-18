#' Takes a sqllite table from an individual caller in maf format
#' and filters as per given user filtering options
#' @param table a maf format dataframe
#' @param impact_values list of IMPACT [optional] values to filter
#' @param hotspot_amino_acid_position_df  [optional] dataframe with Hugo_Sybol, Amino_Acid_Position
#' and hotspot_database columns
#' @param hotspot_genomic_site_df [optional] dataframe with Hugo_Symbol and Chromosome, Start_Position ,
#' End_Position, Hugo_Symbol and hotspot_database columns
#' @return A dataframe with base columns from maf filtered as per user input for Variant_Classification
#' column with "variant_classification" OR IMPACT column with "impact_values" OR Hugo_Symbol using "gene_table".
#' Additional filtering can be done per amino acid position hotspots using "hotspot_amino_acid_site_df"
#' OR per genomic site hotspots using "hotspot_genomic_site_df"
#'

filterMaf <- function(table,
                      impact_values,
                      hotspot_amino_acid_position_df,
                      hotspot_genomic_site_df) {
  if (!missing(impact_values)) {
    # IMPACT filtering values
    table <- table %>%
      dplyr::filter(IMPACT %in% impact_values)
  }

  # filtered calls_base
  calls_base <- table %>%
    as.data.frame() %>%
    dplyr::mutate(
      Amino_Acid_Position =
        # format of Protein_position is 594/766
      # [Amino_Acid_Position/ Protein_length]
      # stripping protein length here to get Amino_Acid_Position
      gsub("/.*", "", Protein_position)
    ) %>%
    mutate(
      Start_Position = as.numeric(Start_Position),
      End_Position = as.numeric(End_Position)
    )

  if (!missing(hotspot_amino_acid_position_df)) {
    # filter by amino acid hotspot_database
    calls_base_aa_filt <- calls_base %>%
      dplyr::inner_join(hotspot_amino_acid_position_df,
        by = c("Amino_Acid_Position", "Hugo_Symbol")
      ) %>%
      dplyr::select(-hotspot_database)
  }

  if (!missing(hotspot_genomic_site_df)) {
    # calls_base genomic ranges
    calls_base_gr <- GenomicRanges::makeGRangesFromDataFrame(calls_base,
      keep.extra.columns = TRUE,
      seqnames.field = "Chromosome",
      start.field = "Start_position",
      end.field = "End_position",
      strand.field = "+"
    )
    # genomic hotspot genomic ranges
    genomic_hotspot_gr <- GenomicRanges::makeGRangesFromDataFrame(hotspot_genomic_site_df,
      keep.extra.columns = TRUE,
      seqnames.field = "Chromosome",
      start.field = "Start_position",
      end.field = "End_position",
      strand.field = "+"
    )

    # subset by genomic region hotspot_database
    calls_base_g <- IRanges::mergeByOverlaps(calls_base_gr, genomic_hotspot_gr)

    calls_base_g_filt <- data.frame(
      Chromosome = as.character(GenomicRanges::seqnames(calls_base_g$calls_base_gr)),
      Start_Position = as.integer(GenomicRanges::start(calls_base_g$calls_base_gr)),
      End_Position = as.integer(GenomicRanges::end(calls_base_g$calls_base_gr))
    ) %>%
      bind_cols(as.data.frame(calls_base_g[, -which(colnames(calls_base_g) %in% c(
        "Chromosome",
        "Start_Position",
        "End_Position",
        "calls_base_gr",
        "genomic_hotspot_gr",
        "hotspot_database"
      ))])) %>%
      dplyr::select(-"Hugo_Symbol.1") %>%
      unique()
  }

  if (!missing(hotspot_genomic_site_df) & !missing(hotspot_amino_acid_position_df)) {
    # if both amino acid and genomic region hotspot are provided
    calls_base_combined <- bind_rows(
      calls_base_aa_filt,
      calls_base_g_filt
    )
    return(calls_base_combined)
  } else if (!missing(hotspot_genomic_site_df)) {
    # if only genomic region hotspot is provided
    return(calls_base_g_filt)
  } else if (!missing(hotspot_amino_acid_position_df)) {
    # if only amino acid hotspot is provided
    return(calls_base_aa_filt)
  } else {
    # if no hotspot list/region is provided return base filtered calls
    return(calls_base)
  }
}
