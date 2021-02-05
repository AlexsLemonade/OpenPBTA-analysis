#' Takes a sqllite table from an individual caller in maf format
#' with a Variant_Classification and gene list to filter.
#' @param con sqllite db connection with strelka, mutect, vardict, lancet caller maf tables
#' @param variant_classification_values list of Variant_Classification values to filter
#' @param impact_values list of IMPACT values to filter
#' @param gene_table gene df with columns "Hugo_Symbol" and "type"
#' @return dataframe with counts per gene and start_position
#' and end_position

getRecurrentStatus<- function(con,
                              variant_classification_values='Missense_Mutation|Splice_Region|In_Frame_Del|Frame_Shift_Del|Splice_Site|Splice_Region|Nonsense_Mutation|Nonstop_Mutation|In_Frame_Ins|Frameshift_Ins',
                              impact_values='HIGH|MODERATE',
                              gene_table){
  # gene list to filter sql lite table
  gene_list <- gene_table$Hugo_Symbol
  # select columns to only use columns needed for a base maf
  columns <-c("Chromosome","Start_Position","End_Position","Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode" ,"Protein_position")
  # Variant_Classification filtering values
  variant_classification_values <- unlist(str_split(variant_classification_values,pattern = "[|]"))
  # IMPACT filtering values
  impact_values <- unlist(str_split(impact_values,pattern = "[|]"))
  
  strelka <- dplyr::tbl(con,"strelka")
  mutect <- dplyr::tbl(con,"mutect")
  vardict <- dplyr::tbl(con,"vardict")
  lancet <- dplyr::tbl(con,"lancet")

  strelka_calls_base <- strelka %>%
      dplyr::filter(Variant_Classification %in% variant_classification_values,
             Hugo_Symbol %in% gene_list) %>%
    dplyr::select(columns) %>%
    as.data.frame()
  
   mutect_calls_base <- mutect %>%
      dplyr::filter(Variant_Classification %in% variant_classification_values,
             Hugo_Symbol %in% gene_list) %>%
   dplyr::select(columns) %>%
   as.data.frame()

  vardict_calls_base <- vardict %>%
      dplyr::filter(Variant_Classification %in% variant_classification_values,
             Hugo_Symbol %in% gene_list) %>%
    dplyr::select(columns) %>%
    as.data.frame()

 lancet_calls_base <- lancet %>%
      dplyr::filter(Variant_Classification %in% variant_classification_values,
             Hugo_Symbol %in% gene_list) %>%
    dplyr::select(columns) %>%
    as.data.frame()

  calls_base <- bind_rows(strelka_calls_base,
			  mutect_calls_base,
			  vardict_calls_base,
			  lancet_calls_base) %>%
  unique()

  calls_base_recurrence <- calls_base %>%
    left_join(gene_table,by="Hugo_Symbol") %>%
    count(Chromosome,Start_Position,End_Position,Protein_position,Hugo_Symbol,type,sort = TRUE)

    
  return(calls_base_recurrence)
  
}

