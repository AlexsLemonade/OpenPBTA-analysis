#' @description CNS tumors have subtypes as per the [WHO 2016 CNS subtypes](https://link.springer.com/content/pdf/10.1007/s00401-016-1545-1.pdf). 
#' However, these are not captured in our molecular data so would need to be updated by
#' searching for terms in the reported pathology_free_text_diagnosis column in OpenPBTA histology file.
#' recode_integrated_diagnosis() can be used to recode a integreated_diagnosis based on the
#' pathology_free_text_diagnosis in addition to molecular subtypes.
#' @param histologies_df : Dataframe with Kids_First_Biospecimen_ID,Kids_First_Participant_ID,
#' sample_id,tumor_descriptor,pathology_free_text_diagnosis,broad_histology,short_histology,
#' molecular_subtype,integrated_diagnosis
#' @param pathology_diagnois : A term to filter histologies_df to samples that will be recoded
#' @param include_path_free_text_dx_terms : Regex term to match with values 
#' in pathology_free_text_diagnosis in histologies_df to be included in the recoding analysis
#' @param exclude_path_free_text_dx_terms : Default = NULL; Regex term to match with values in 
#' pathology_free_text_diagnosis in histologies_df to be excluded from the recoding analysis
#' @param old_integrated_diagnosis_term : Old term in integrated_diagnosis that should be replaced 
#' by replace_integrated_diagnosis_term
#' @param replace_integrated_diagnosis_term : New term for integrated_diagnosis that will replace 
#' an old integrated_term


recode_integrated_diagnosis <- function(histologies_df,
                                        pathology_diagnosis,
                                        include_path_free_text_dx_terms,
                                        exclude_path_free_text_dx_terms=NULL,
                                        old_integrated_diagnosis_term,
                                        replace_integrated_diagnosis_term){
  subtype_include_df <- histologies_df %>%
    # We are only concerned with samples where the pathology_free_text_diagnosis
    # contains the include term
    dplyr::filter(pathology_diagnosis == pathology_diagnosis,
                  grepl(include_path_free_text_dx_terms,pathology_free_text_diagnosis)) %>%
    # Retain only relevant identifier and disease label columns
    dplyr::select(Kids_First_Biospecimen_ID,
                  Kids_First_Participant_ID,
                  sample_id,
                  tumor_descriptor,
                  pathology_free_text_diagnosis,
                  broad_histology,
                  short_histology,
                  molecular_subtype,
                  integrated_diagnosis) %>%
    # To smooth the way for string detection for the pathology free text, we
    # add a column where all of the text is converted to lowercase
    dplyr::mutate(pathology_free_text_dx_lower = 
                    stringr::str_to_lower(pathology_free_text_diagnosis)) %>%
    # String detection in pathology free text per the table
    dplyr::mutate(
      integrated_diagnosis = dplyr::case_when(
        # if NA integrated diagnosis keep as NA since there is no subtype to integrate
        is.na(integrated_diagnosis) ~ NA_character_,
        # if there is an existing integrated_diagnosis then replace with new integrated_diagnosis
        stringr::str_detect(pathology_free_text_dx_lower,
                            include_path_free_text_dx_terms) ~ str_replace(integrated_diagnosis,old_integrated_diagnosis_term,replace_integrated_diagnosis_term),
      ),
      harmonized_diagnosis = dplyr::case_when(
        # if NA integrated diagnosis keep as NA since there is no subtype to integrate
        is.na(integrated_diagnosis) ~ replace_integrated_diagnosis_term,
        # if there is an existing integrated_diagnosis then replace with new integrated_diagnosis
        stringr::str_detect(pathology_free_text_dx_lower,
                            include_path_free_text_dx_terms) ~ str_replace(integrated_diagnosis,old_integrated_diagnosis_term,replace_integrated_diagnosis_term),
      ),
      Notes = if_else(is.na(integrated_diagnosis),"Updated via pathology_free_text_diagnosis" , "Updated via OpenPBTA subtyping and pathology_free_text_diagnosis")
    ) 
  
  if(!is.null(exclude_path_free_text_dx_terms)){
    subtype_include_df <- subtype_include_df %>%
      filter(!grepl(exclude_path_free_text_dx_terms,pathology_free_text_diagnosis))
    
  }
  
  get_subtype_df<- subtype_include_df %>%
    # Drop the column we added for convenience of string detection
    # and to format to match compiled_mol_subtypes_pathology_clinical_df
    dplyr::select(-pathology_free_text_dx_lower,
                  -pathology_free_text_diagnosis)
  
  return(get_subtype_df)
}
