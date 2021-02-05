#' Takes a maf filename and adds hotspot status from 
#' MSKCC snv/indel cancer hotspot, and internal hotspot
#' dataframes with required column "Hugo_Symbol" and 
#' "Amino_Acid_Postion".
#' @param maf_file file name should have "Chromosome","Start_Position",
#' "End_Position","Hugo_Symbol","Protein_position" columns
#' @return dataframe with hotspot status annotated as "Yes" 
#' if found in the given cancer and internal hotspot dataframes

getHotspotStatus<- function(maf_file){
  calls_base_hotspot <- data.table::fread(maf_file ,
                                          select = c("Chromosome",
                                                     "Start_Position",
                                                     "End_Position",
                                                     "Hugo_Symbol",
                                                     "Protein_position"),
                                          data.table = FALSE) %>%
    unique() %>%
    dplyr::mutate(
      Amino_Acid_Position = 
        # strip protein length
        gsub("/.*","",Protein_position)
    ) %>%
    dplyr::mutate(
      hotspot = case_when(
        Amino_Acid_Position %in% 
          # if overlaps the hotspot snv Amino_Acid_Position 
          hotspot_database_2017_snv$Amino_Acid_Position &
          Hugo_Symbol %in% hotspot_database_2017_snv$Hugo_Symbol ~ "Yes",
        Amino_Acid_Position %in% 
          # if overlaps the hotspot snv Amino_Acid_Position 
          hotspot_database_2017_indel$Amino_Acid_Position &
          Hugo_Symbol %in%  hotspot_database_2017_indel$Hugo_Symbol~ "Yes",
        Amino_Acid_Position %in% 
          # if overlaps the hotspot snv Amino_Acid_Position 
          internal_hotspot$Amino_Acid_Position &
          Hugo_Symbol %in% internal_hotspot$Hugo_Symbol~ "Yes"
      )
    ) %>%
    dplyr::filter(hotspot == "Yes")
  
  return(calls_base_hotspot)
  
}
