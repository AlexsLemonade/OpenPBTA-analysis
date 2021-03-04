#' Takes a sqllite table from an individual caller in maf format
#' and filters as per given user filtering options
#' @param table prefered sqllite table with individual caller maf or a maf format dataframe
#' @param variant_classification_values [optional] list of Variant_Classification values to filter
#' @param impact_values list of IMPACT [optional] values to filter
#' @param additional_cols [optional] character vector of maf columns to add other than c("Chromosome",
#'                                                                            "Start_Position",
#'                                                                            "End_Position",
#'                                                                            "Reference_Allele",
#'                                                                            "Hugo_Symbol",
#'                                                                            "Variant_Classification",
#'                                                                            "IMPACT",
#'                                                                            "Tumor_Sample_Barcode" ,
#'                                                                            "Protein_position)
#' @param gene_table [optional] gene df with columns "Hugo_Symbol" and "type"
#' @param hotspot_amino_acid_position_df  [optional] dataframe with Hugo_Sybol, Amino_Acid_Position and hotspot_database columns
#' @param hotspot_genomic_site_df [optional] dataframe with Hugo_Symbol and Chromosome, Start_Position , 
#' End_Position and hotspot_database columns 
#' @return A dataframe with base columns from maf filtered as per user input for Variant_Classification 
#' column with "variant_classification" OR IMPACT column with "impact_values" OR Hugo_Symbol using "gene_table".
#' Additional filtering can be done per amino acid position hotspots using "hotspot_amino_acid_site_df" 
#' OR per genomic site hotspots using "hotspot_genomic_site_df"
#'

prepMaf <- function(table,
                    variant_classification_values,
                    impact_values,
                    gene_table,
                    additional_cols = "VAF,HGVSp_Short,gnomAD_AF,Variant_Type,Tumor_Seq_Allele2,dbSNP_RS",
                    hotspot_amino_acid_position_df, 
                    hotspot_genomic_site_df  ) {
  if (!missing(gene_table)){
  # gene list to filter sql lite table
  gene_list <- gene_table$Hugo_Symbol
  table <- table %>%
    dplyr::filter(Hugo_Symbol %in% gene_list)
  }
  
  # base maf columns
  columns <- c("Chromosome",
                     "Start_Position",
                     "End_Position",
                     "Reference_Allele",
                     "Hugo_Symbol",
                     "Variant_Classification",
                     "IMPACT",
                     "Tumor_Sample_Barcode" ,
                     "Protein_position")
  
 
  # select columns to only use columns needed for a base maf
  additional_cols = unlist(str_split(additional_cols,pattern = "[,]"))
  columns <- unique(c( columns ,
                     additional_cols))
  
  
  if (!missing(variant_classification_values)){
  # Variant_Classification filtering values
  variant_classification_values <- unlist(str_split(variant_classification_values,pattern = "[|]"))
  table <- table %>%
    dplyr::filter(Variant_Classification %in% variant_classification_values)
  }
  
  if (!missing(impact_values)){
  # IMPACT filtering values
  impact_values <- unlist(str_split(impact_values,pattern = "[|]"))
  table <- table %>%
    dplyr::filter(IMPACT %in% impact_values)
  }
  
  # filtered calls_base 
  calls_base <- table %>%
    dplyr::select(!!!(columns)) %>%
    as.data.frame() %>%
    dplyr::mutate(
      Amino_Acid_Position =
        # strip protein length
        gsub("/.*","",Protein_position)
    )
  
  if (!missing(hotspot_amino_acid_position_df)){
    # filter by amino acid hotspot_database
    calls_base_aa_filt <- calls_base %>%
    left_join(hotspot_amino_acid_position_df,by=c("Amino_Acid_Position","Hugo_Symbol")) %>%
      dplyr::filter(!is.na(hotspot_database)) %>%
      dplyr::select(-hotspot_database)
  } 
  
  if(!missing(hotspot_genomic_site_df)){
    # calls_base genomic ranges
    calls_base_gr <- GenomicRanges::makeGRangesFromDataFrame(calls_base,keep.extra.columns = TRUE,
                                                             seqnames.field = "Chromosome",
                                                             start.field = "Start_position",
                                                             end.field = "End_position",
                                                             strand.field = "+")
    # genomic hotspot genomic ranges
    genomic_hotspot_gr <- GenomicRanges::makeGRangesFromDataFrame(hotspot_genomic_site_df,keep.extra.columns = TRUE,
                                                                  seqnames.field = "Chromosome",
                                                                  start.field = "Start_position",
                                                                  end.field = "End_position",
                                                                  strand.field = "+") 
    
    # subset by genomic region hotspot_database
    calls_base_g <- IRanges::mergeByOverlaps(calls_base_gr,genomic_hotspot_gr)
    calls_base_g_filt <- data.frame(Chromosome= as.character( GenomicRanges::seqnames(calls_base_g$calls_base_gr)),
                             Start_Position=as.integer( GenomicRanges::start(calls_base_g$calls_base_gr)),
                             End_Position=as.integer( GenomicRanges::end(calls_base_g$calls_base_gr)),
                             Hugo_Symbol = calls_base_g$Hugo_Symbol,
                             Variant_Classification = calls_base_g$Variant_Classification,
                             IMPACT = calls_base_g$IMPACT,
                             Protein_position = calls_base_g$Protein_position,
                             Amino_Acid_Position = calls_base_g$Amino_Acid_Position,
                             Reference_Allele= calls_base_g$Reference_Allele,
                             Variant_Type = calls_base_g$Variant_Type,
                             Tumor_Seq_Allele2 = calls_base_g$Tumor_Seq_Allele2,
                             Tumor_Sample_Barcode = calls_base_g$Tumor_Sample_Barcode,
                             HGVSp_Short = calls_base_g$HGVSp_Short,
                             VAF = calls_base_g$VAF) %>%
      unique() 
  }
  
  if(!missing(hotspot_genomic_site_df) & !missing(hotspot_amino_acid_position_df)){
    # if both amino acid and genomic region hotspot are provided
    calls_base_combined <- bind_rows(calls_base_aa_filt ,
                                     calls_base_g_filt)
    return(calls_base_combined)
  } else if(!missing(hotspot_genomic_site_df) & missing(hotspot_amino_acid_position_df)){
    # if only genomic region hotspot is provided
    return(calls_base_g_filt)
  } else if(missing(hotspot_genomic_site_df) & !missing(hotspot_amino_acid_position_df)){
    # if only amino acid hotspot is provided
    return(calls_base_aa_filt)
  } else {
    # if no hotspot list/region is provided return base filtered calls
    return(calls_base)
  }
  
  
}


