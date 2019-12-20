# independent_samples.R

#' Generate a vector of unique samples
#' 
#' The samples from this function will be unique with respect to participants
#' i.e. only no two samples will come from the same participant. The input list 
#' should be pre-filtered by `experimental_strategy` and `sample_type`.  
#' 
#' 
#' @param sample_df A data frame of samples, with columns corresponding to those
#'   in `pbta-histologies.tsv`
#' @param tumor_types Designates which types of tumors will be included. Options
#'   are "primary" to include only primary tumors, "prefer_primary" to include
#'   primary tumors when available, but fall back to other types, or "any" to
#'   randomly select among all available specimens. As of v5, primary tumors
#'   are defined as those designated "Initial CNS Tumor" or "Diagnosis" in the
#'   `tumor_descriptor` field.
#' @param seed An optional random number seed. 
#' 
#' @return a data frame of Participant and Specimen IDs, each present only once.
independent_samples <- function(sample_df, 
                                tumor_types = c("primary", "prefer_primary", "any"), 
                                seed){
  tumor_types <- match.arg(tumor_types)
  if(!missing(seed)){set.seed(seed)}
  
  primary_descs <- c("Initial CNS Tumor", "Diagnosis")
  
  if(tumor_types == "prefer_primary"){
    # find cases where non-primary is the only option
    no_primary <- sample_df %>% 
      dplyr::group_by(Kids_First_Participant_ID) %>%
      dplyr::summarize(n_primary = sum(tumor_descriptor %in% primary_descs)) %>%
      dplyr::filter(n_primary == 0) %>%
      dplyr::pull(Kids_First_Participant_ID)
  } else {
    no_primary <- c()
  }
  
  if(tumor_types %in% c("primary", "prefer_primary")){
    primary_df <- sample_df %>%
      dplyr::filter(tumor_descriptor %in% primary_descs)
    noprimary_df <- sample_df %>%
      dplyr::filter(Kids_First_Participant_ID %in% no_primary)
    sample_df <- dplyr::bind_rows(primary_df, noprimary_df)
  } 

  # get the samples from the earliest timepoints for each Participant
  # age_at_diagnosis_days is no longer relevant,
  # as it is the same for all samples from an participant, but
  # leaving this in for future use in case we get specimen order data 
  early_samples <- sample_df %>%
   dplyr::group_by(Kids_First_Participant_ID) %>%
   dplyr::summarize(age_at_diagnosis_days = min(age_at_diagnosis_days)) %>%
   dplyr::left_join(sample_df, by = c("Kids_First_Participant_ID",
                                      "age_at_diagnosis_days"))

  # Choose randomly among specimens from the same participant
  early_ind <- early_samples %>%
    dplyr::group_by(Kids_First_Participant_ID) %>%
    dplyr::summarize(Kids_First_Biospecimen_ID = sample(Kids_First_Biospecimen_ID, 1)) 
  
  return(early_ind)
}
