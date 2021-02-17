#' Takes a sqllite table from an individual caller in maf format
#' with a Variant_Classification and gene list to filter.
#' @param table prefered sqllite table with individual caller maf or a maf format dataframe
#' @param variant_classification_values list of Variant_Classification values to filter
#' @param impact_values list of IMPACT values to filter
#' @param gene_table gene df with columns "Hugo_Symbol" and "type"
#' @return dataframe with base columns from maf
#'
#'

prepMaf <- function(table,
                              variant_classification_values='Missense_Mutation|Splice_Region|In_Frame_Del|Frame_Shift_Del|Splice_Site|Splice_Region|Nonsense_Mutation|Nonstop_Mutation|In_Frame_Ins|Frameshift_Ins',
                              impact_values='HIGH|MODERATE',
                              gene_table,
                    additional_cols = "VAF,HGVSp_Short,gnomAD_AF,Variant_Type,Tumor_Seq_Allele2,dbSNP_RS"){
  # gene list to filter sql lite table
  gene_list <- gene_table$Hugo_Symbol
  # select columns to only use columns needed for a base maf
  additional_cols = unlist(str_split(additional_cols,pattern = "[,]"))
  columns <-unique(c("Chromosome",
                     "Start_Position",
                     "End_Position",
                     "Reference_Allele",
                     "Hugo_Symbol",
                     "Variant_Classification",
                     "IMPACT",
                     "Tumor_Sample_Barcode" ,
                     "Protein_position", 
                     additional_cols))
  # Variant_Classification filtering values
  variant_classification_values <- unlist(str_split(variant_classification_values,pattern = "[|]"))
  # IMPACT filtering values
  impact_values <- unlist(str_split(impact_values,pattern = "[|]"))
  
  
  calls_base <- table %>%
      dplyr::filter(Variant_Classification %in% variant_classification_values,
             Hugo_Symbol %in% gene_list,
             IMPACT %in% impact_values) %>%
    dplyr::select(columns) %>%
    as.data.frame() %>%
    dplyr::mutate(
      Amino_Acid_Position =
        # strip protein length
        gsub("/.*","",Protein_position)
    )

  return(calls_base)
  
}


