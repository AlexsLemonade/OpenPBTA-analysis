# independent_samples.R

#' Generate a vector of unique samples
#' 
#' The samples from this function will be unique with respect to participants
#' i.e. only no two samples will come from the same participant. The input list 
#' should be pre-filtered by `experimental_strategy` and `sample_type`.  
#' 
#' 
#' @param histology_df A data frame of samples, with columns corresponding to those
#'   in `histologies.tsv`
#' @param independent_level Designates whether we want to count independent samples in 
#'  different cohorts as independent or not. "all-cohorts" consider the same sampe
#'  in different cohorts as the same sample and "each-cohort" consider the same sample
#'  in different cohorts as "independent" (different). 
#' @param tumor_types Designates which types of tumors will be included. Options
#'   are "primary" to include only primary tumors, "prefer_primary" to include
#'   primary tumors when available, but fall back to other types, or "any" to
#'   randomly select among all available specimens. As of v6, primary tumors
#'   are defined as those designated "Initial CNS Tumor" or "Primary Tumor" in the
#'   `tumor_descriptor` field.
#' @param seed An optional random number seed. 
#' 
#' @return a data frame of Participant and Specimen IDs, each present only once.
independent_samples <- function(histology_df, 
                                tumor_types = c("primary", "relapse", "prefer_primary", "any"), 
                                independent_level = c("all-cohorts", "each-cohort"),
                                seed){
  tumor_types <- match.arg(tumor_types)
  independent_level <- match.arg(independent_level)
  if(!missing(seed)){set.seed(seed)}
  
  primary_descs <- c("Initial CNS Tumor", "Primary Tumor")
  relapse_descs <- c("Recurrence", "Progressive", "Progressive Disease Post Mortem")
  
  if(tumor_types %in% c("prefer_primary")){
    # find cases where non-primary is the only option
    no_primary <- histology_df %>% 
      dplyr::group_by(Kids_First_Participant_ID) %>%
      dplyr::summarize(n_primary = sum(tumor_descriptor %in% primary_descs)) %>%
      dplyr::filter(n_primary == 0) %>%
      dplyr::pull(Kids_First_Participant_ID)
  } else {
    no_primary <- c()
  }
  
  if(tumor_types %in% c("primary", "prefer_primary")){
    primary_df <- histology_df %>%
      dplyr::filter(tumor_descriptor %in% primary_descs)
    noprimary_df <- histology_df %>%
      dplyr::filter(Kids_First_Participant_ID %in% no_primary)
    sample_df <- dplyr::bind_rows(primary_df, noprimary_df)
  } 
  
  if(tumor_types == "relapse"){
    sample_df <- histology_df %>%
      dplyr::filter(tumor_descriptor %in% relapse_descs)
  } 

  if(independent_level == "each-cohort"){
    # find out the list of cohorts and cancer groups and loop through them to run the process for each cancer group and cohort
    cohort_list <- sample_df$cohort %>% unique()
    cancer_group_list <- sample_df$cancer_group %>% unique()   
    
    independent_each <- data.frame(Kids_First_Participant_ID = character(), Kids_First_Biospecimen_ID = character(), stringsAsFactors = FALSE)
   
    # loop through each cohort and cancer group and combine the results together
    for (i in 1:length(cohort_list)){
      for (j in 1:length(cancer_group_list)){
        # filter to the specific cancer group and cohort
        cohort_name <- cohort_list[i]
        cancer_group_name <- cancer_group_list[j]
        # deal with cancer group is NA to avoid missing samples
        if(is.na(cancer_group_name)){
          filtered_df <- sample_df %>% dplyr::filter(cohort == cohort_name) %>%
            dplyr::filter(is.na(cancer_group))
        }else{
          filtered_df <- sample_df %>% dplyr::filter(cohort == cohort_name) %>%
            dplyr::filter(cancer_group == cancer_group_name)
          }
        # some specific group does not have any specimen 
        if(nrow(filtered_df) == 0){
          independent_filtered <- data.frame(Kids_First_Participant_ID = character(), Kids_First_Biospecimen_ID = character(), stringsAsFactors = FALSE)
          independent_each <- rbind(independent_each, independent_filtered)
        }else{
          # find the independent samples for the specific cancer group and cohort
          independent_filtered <- filtered_df %>%
            dplyr::group_by(Kids_First_Participant_ID) %>%
            dplyr::summarize(Kids_First_Biospecimen_ID = sample(Kids_First_Biospecimen_ID, 1)) %>%
            data.frame() %>% dplyr::distinct()
          # merge the independent samples together
          independent_each <- rbind(independent_each, independent_filtered)
        }
      }
    }
    return(independent_each)
  }
  
  
  if(independent_level == "all-cohorts"){
    independent_all<- sample_df %>%
      dplyr::group_by(Kids_First_Participant_ID) %>%
      dplyr::summarize(Kids_First_Biospecimen_ID = sample(Kids_First_Biospecimen_ID, 1)) %>%
      data.frame() %>% dplyr::distinct()
    return(independent_all)
  }
}
