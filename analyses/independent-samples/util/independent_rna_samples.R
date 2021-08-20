# independent_rna_samples.R

#' Generate a vector of unique rna samples
#' 
#' The samples from this function will be unique with respect to participants
#' i.e. only no two samples will come from the same participant. The input list 
#' should be pre-filtered by `composition` and `sample_type`.  
#' 
#' 
#' @param independent_dna_sample_df A data frame of samples, with columns 
#' corresponding to those in `independent-specimens.wgswxspanel.primary.tsv` 
#' or `independent-specimens.wgswxspanel.relapse.tsv` 
#' or `independent-specimens.wgswxspanel.primary-plus.tsv` depending on what 
#' set of samples you need to include
#' @param histology_df A data frame of samples, with columns corresponding 
#' to those `histologies.tsv`
#' @param independent_level Designates whether we want to count independent samples in 
#' different cohorts as independent or not. "all-cohorts" consider the same sampe
#' in different cohorts as the same sample and "each-cohort" consider the same sample
#' in different cohorts as "independent" (different). 
#' @param match_type Designates which type matching needs to be done. Options
#' are "independent_dna" to include only rna samples that match the 
#' independent-specimens.wgswxspanel sample set, 
#' "independent_dna_plus_only_rna" to include samples that macth the dna sample 
#' set plus include samples where only rna samples exists
#' @param tumor_description_rna_only Tumor descriptors to select samples where
#' only RNA samples are available and will have no matching id in independent_dna_sample_df
#' Opetions are "primary" to select only primary/initial tumors. Primary tumors are defined as those designated "Initial CNS Tumor"+ "Primary Tumor".
#' "primary_plus" if you would like to select other non-initial tumor RNA-Seq sample if no 
#' initial tumor RNA-Seq sample exists
#' or "relapse" in select only relapse tumors. Relapse tumors are defined as those designated by "Recurrence", "Recurrent Tumor","Recurrent tumor","Progressive","Progressive Disease Post-Mortem" in `tumor_descriptor` field
#' @param seed An optional random number seed. 
#' 
#' @return a data frame of Participant and Specimen IDs, each present only once.
independent_rna_samples <- function(independent_dna_sample_df, 
                                histology_df,
                                independent_level = c("all-cohorts", "each-cohort"),
                                match_type = c("independent_dna", "independent_dna_plus_only_rna"),
                                tumor_description_rna_only = c("primary", "relapse", "primary_plus"),
                                seed){
  match_type <- match.arg(match_type)
  tumor_description_rna_only <- match.arg(tumor_description_rna_only)
  independent_level <- match.arg(independent_level)
  if(!missing(seed)){set.seed(seed)}
  primary_descs <- c("Initial CNS Tumor", "Primary Tumor")
  relapse_descs <- c("Recurrence", "Progressive", "Progressive Disease Post Mortem")
  
  # Find sample set for the dna independent samples 
  # This will always be the included since in both the following
  # conditions "independent_dna" "independent_dna_plus_only_rna"
  # 
  independent_dna <- histology_df %>%
    # include matched independent_dna samples
    dplyr::filter(Kids_First_Participant_ID %in%
                    independent_dna_sample_df$Kids_First_Participant_ID)
  matched_rna <- histology_df %>%
    # keep rna from histology_df
    dplyr::filter(experimental_strategy == "RNA-Seq",
                  # find participants which have matching dna samples in independent_wgswxspanel
                  Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
                  # keep specific sample_ids since some participants might have multiple sample_ids
                  sample_id %in% independent_dna$sample_id)
  
  matched_rna_primary <- histology_df %>%
    # keep rna from histology_df
    dplyr::filter(experimental_strategy == "RNA-Seq",
                  tumor_descriptor %in% primary_descs,
                  Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
                  sample_id %in% independent_dna$sample_id)
  
  matched_rna_relapse <- histology_df %>%
    # keep rna from histology_df
    dplyr::filter(experimental_strategy == "RNA-Seq",
                  tumor_descriptor %in% relapse_descs,
                  Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
                  sample_id %in% independent_dna$sample_id)

  # Here we are adding only initial only-RNA-Seq samples
  # since this will always to part of independent_dna_plus_only_rna
  # regardless tumor_description_rna_only is "primary" OR "primary_plus"
  #
  if( match_type == "independent_dna_plus_only_rna" & tumor_description_rna_only == "primary") {
    # find sample set where we initial only-RNA-Seq samples
    only_rna_initial <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq",
                    tumor_descriptor %in% primary_descs,
                    # find and remove participants which have 
                    # matching dna samples in independent_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
    
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors
    sample_df <- bind_rows(matched_rna_primary, only_rna_initial)
  }
  
  if( match_type == "independent_dna_plus_only_rna" & tumor_description_rna_only == "relapse") {
    # find sample set where we initial only-RNA-Seq samples
    only_rna_relapse <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq",
                    tumor_descriptor %in% relapse_descs,
                    # find and remove participants which have 
                    # matching dna samples in independent_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
    
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors
    sample_df <- bind_rows(matched_rna_relapse,only_rna_relapse)
  }
  
  # Here we are adding only-RNA-Seq samples which are not initial
  # if tumor_description_rna_only == "primary_plus"
  #
  if(match_type == "independent_dna_plus_only_rna" & tumor_description_rna_only == "primary_plus"){
    # find sample set where we only find rna samples
    only_rna_initial <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq",
                    tumor_descriptor %in% primary_descs,
                    # find and remove participants which have 
                    # matching dna samples in independent_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID)
    only_rna_plus <- histology_df %>% 
      # keep rna from histology_df
      dplyr::filter(experimental_strategy == "RNA-Seq",
                    # find and remove participants which have 
                    # matching dna samples in independend_wgswxspanel
                    !Kids_First_Participant_ID %in% independent_dna$Kids_First_Participant_ID,
                    # and participant not in only_rna_initial sample set
                    !Kids_First_Participant_ID %in% only_rna_initial$Kids_First_Participant_ID
                    )
    # has rna samples which match the independent samples provided plus rna only sample which are primary tumors plus rna samples where no primary primaries exists
    sample_df <- bind_rows(matched_rna,only_rna_initial, only_rna_plus)
  } 
  
  if(independent_level == "each-cohort"){
    # find out the list of cohorts and cancer groups and loop through them to run the process for each cancer group and cohort
    # combination of cohort and cancer groups
    cohort_cancer_group_combo <- sample_df %>%
      dplyr::select(cohort, cancer_group) %>%
      unique() 
    
    independent_each <- data.frame(Kids_First_Participant_ID = character(), 
                                   Kids_First_Biospecimen_ID = character(), 
                                   stringsAsFactors = FALSE)
    
    # loop through each cohort and cancer group and combine the results together
    for(i in 1:nrow(cohort_cancer_group_combo)){
      # filter to the specific cancer group and cohort
      cohort_name <- cohort_cancer_group_combo[i,1] %>% as.character()
      cancer_group_name <- cohort_cancer_group_combo[i,2] %>% as.character()
      
      # deal with cancer group is NA to avoid missing samples
      if(is.na(cancer_group_name)){
        filtered_df <- sample_df %>% 
          dplyr::filter(cohort == cohort_name) %>%
          dplyr::filter(is.na(cancer_group))
      }else{
        filtered_df <- sample_df %>% 
          dplyr::filter(cohort == cohort_name) %>%
          dplyr::filter(cancer_group == cancer_group_name)
      }
      
      # find the independent samples for the specific cancer group and cohort
      # "If there are multiple rows for a given combination of inputs, only the first
      # row will be preserved. If omitted, will use all variables." -- distinct in dplyr 0.8.3
      independent_filtered <- filtered_df %>%
        dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>%
        dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)
      
      # merge the independent samples together
      independent_each <- rbind(independent_each, independent_filtered)
    }
    return(independent_each)
    }
  
  if(independent_level == "all-cohorts"){
    
    independent_all <- sample_df %>%
      dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>%
      dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID)
    
    return(independent_all)
  }
}

