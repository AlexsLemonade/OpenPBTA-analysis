# independent_rna_samples.R

#' Generate a vector of unique rna samples
#' 
#' The samples from this function will be unique with respect to participants
#' i.e. only no two samples will come from the same participant. The input list 
#' should be pre-filtered by `composition` and `sample_type`.  
#' 
#' 
#' @param independent_dna_sample_df A data frame of samples, with columns 
#' corresponding to those in `independent-specimens.wgswxs.primary.tsv` 
#' or `independent-specimens.wgswxs.primary-plus.tsv` depending on what 
#' set of samples you need to include
#' @param histology_df A data frame of samples, with columns corresponding 
#' to those `pbta-histologies.tsv`
#' @param match_type Designates which type matching needs to be done. Options
#' are "independent_dna" to include only rna samples that match the 
#' independent-specimens.wgswxs sample set, 
#' "independent_dna_plus_only_rna" to include samples that macth the dna sample 
#' set plus include samples where only rna samples exists
#' @param tumor_description_rna_only Tumor descriptors to select samples where
#' only RNA samples are available and will have no matching id in independent_dna_sample_df
#' Opetions are "primary" to select only primary/initial tumors As of v18, primary tumors are defined as those designated "Initial CNS Tumor" .
#' "primary_plus" if you would like to select other non-initial tumor RNA-Seq sample if no 
#' initial tumor RNA-Seq sample exists
#' or "Diagnosis" in the `tumor_descriptor` field.
#' @param seed An optional random number seed. 
#' 
#' @return a data frame of Participant and Specimen IDs, each present only once.
independent_rna_samples <- function(independent_dna_sample_df, 
                                histology_df,
                                match_type = c("independent_dna", "independent_dna_plus_only_rna"),
                                tumor_description_rna_only = c("primary","primary_plus"),
                                seed){
  match_type <- match.arg(match_type)
  tumor_description_rna_only <- match.arg(tumor_description_rna_only)
  if(!missing(seed)){set.seed(seed)}
  primary_descs <- c("Initial CNS Tumor", "Diagnosis")
  
  
  if(match_type == "independent_dna" | match_type == "independent_dna_plus_only_rna"){
    # find sample set for the dna independent samples 
    independent_dna <- histology_df %>%
      dplyr::filter(Kids_First_Biospecimen_ID %in%
                      independent_dna_sample_df$Kids_First_Biospecimen_ID)
    matched_rna <- histology_df %>%
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq",
                    # find participants which have matching dna samples in independent_wgswxs
                    Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
                    # keep specific sample_ids since some participants might have multiple sample_ids
                    sample_id %in% independent_dna$sample_id)
    
    # has rna samples which match the independent samples provided
    sample_df <- matched_rna
  } 
  
  if(match_type == "independent_dna_plus_only_rna" & 
     (tumor_description_rna_only == "primary" | tumor_description_rna_only == "primary_plus")) {
    # find sample set where we only find rna samples
    only_rna_intial <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq",
                    tumor_descriptor %in% primary_descs,
                    # find and remove participants which have 
                    # matching dna samples in independend_wgswxs 
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
    
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors
    sample_df <- bind_rows(sample_df,only_rna_intial)
  }
  
  if(match_type == "independent_dna_plus_only_rna" & tumor_description_rna_only == "primary_plus"){
    # find sample set where we only find rna samples
    only_rna_plus <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq",
                    # find and remove participants which have 
                    # matching dna samples in independend_wgswxs 
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
                    # and participant not in only_rna_initial sample set
                    !Kids_First_Participant_ID %in% only_rna_intial$Kids_First_Participant_ID
                    )
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors plus rna samples where no primary primaries exists
    sample_df <- bind_rows(sample_df,only_rna_plus)
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
