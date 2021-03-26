#' Takes a table from an individual caller in maf format
#' and filters as per given user filtering options
#' @param table a maf format dataframe
#' @param impact_values vector of IMPACT [optional] values to filter
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
                      hotspot_database_2017_snv_df,
                      hotspot_database_2017_indel_df,
                      hotspot_genomic_site_df) {
  if (!missing(impact_values)) {
    # IMPACT filtering values
    table <- table %>%
      dplyr::filter(IMPACT %in% impact_values)
  }
  
  # filtered calls_base
  calls_base <- table %>%
    # only keep annotation for canonical protein_coding transcripts 
    filter( CANONICAL=='YES' & 
              BIOTYPE=='protein_coding') %>%
    as.data.frame() %>%
    # create an Amino_Acid_Position column
    dplyr::mutate(
      Amino_Acid_Position =
        case_when(
          # If not splice site we want to capture Amino_Acid_Position to match
          !grepl("Splice_Site",Variant_Classification) ~
            # format of Protein_position is 594/766
            # [Amino_Acid_Position/ Protein_length]
            # stripping protein length here to get Amino_Acid_Position
            gsub("/.*", "", Protein_position),
          # If splice site we want to capture HGVSp_Short to match 
          grepl("Splice_Site",Variant_Classification) ~
            # since these sites are intronic sites adding
            # HGSVp_Short instead of Protein_position
            gsub("p[.]","",HGVSp_Short))
    ) %>%
    mutate(
      Start_Position = as.numeric(Start_Position),
      End_Position = as.numeric(End_Position)
    )

  if (!missing(hotspot_database_2017_snv_df) & !missing(hotspot_database_2017_indel_df)) {
    # filter by amino acid hotspot_database
    calls_base_aa_filt_match <- calls_base %>%
      dplyr::inner_join(hotspot_database_2017_snv_df,
        by = c("Amino_Acid_Position", "Hugo_Symbol")
      ) %>%
      select(-Amino_Acid_Position)
    
    # if indel then check for overlap 
    # within region of Amino_Acid_Position
    calls_base_aa_filt_indel <- calls_base %>%
      # if indel in Variant_Classification
      filter(grepl("Ins|Del",Variant_Classification )) %>%
      mutate( Amino_Acid_Start = gsub("-.*" , "", Amino_Acid_Position),
              Amino_Acid_End = gsub ("*.-" , "" , Amino_Acid_Position)
              ) %>%
      dplyr:: left_join(hotspot_database_2017_indel_df,
                        by = c( "Hugo_Symbol"),
                        suffix = c(".calls_base",".hotspot_indel_df")  
      ) %>% 
      # keep only sites that where 
      # calls$start >= hotspot_indel$start and
      # calls$end <= hotspot_indel$stop
      # start and end could be “?” (question mark) which is used to indicate unknown positions as well
      filter(Amino_Acid_Start.calls_base <= Amino_Acid_End.hotspot_indel_df | Amino_Acid_Start.calls_base=="?",
             Amino_Acid_End.calls_base >= Amino_Acid_Start.hotspot_indel_df | Amino_Acid_End.calls_base == "?"
             ) %>%
      select(-starts_with("Amino_Acid_"))
    
    calls_base_aa_filt <- bind_rows( calls_base_aa_filt_match,
                                     calls_base_aa_filt_indel)
    
    
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
    calls_base_g <- IRanges::subsetByOverlaps(calls_base_gr, genomic_hotspot_gr)

    calls_base_g_filt <- as.data.frame(calls_base_g) %>%
      rename(Chromosome=seqnames,
             Start_Position=start,
             End_Position=end) %>%
      select(-width,-strand,-Amino_Acid_Position) %>%
      unique()
  }

  if (!missing(hotspot_genomic_site_df) & !missing(hotspot_database_2017_snv_df) & !missing(hotspot_database_2017_indel_df)) {
    # if both amino acid and genomic region hotspot are provided
    calls_base_combined <- bind_rows(
      calls_base_aa_filt,
      calls_base_g_filt
    )
    return(calls_base_combined)
  } else if (!missing(hotspot_genomic_site_df)) {
    # if only genomic region hotspot is provided
    return(calls_base_g_filt)
  } else if (!missing(hotspot_database_2017_snv_df) & !missing(hotspot_database_2017_indel_df)) {
    # if only amino acid hotspot is provided
    return(calls_base_aa_filt)
  } else {
    # if no hotspot list/region is provided return base filtered calls
    return(calls_base)
  }
}
