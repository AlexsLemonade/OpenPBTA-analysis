# Functions for calculating quantiles of variables used for survival analysis

# using 0.5 as a cutoff
calc_groups <- function(meta_indep, tp53_quantiles, telomerase_quantiles){
  
  # extract the quantile infomration 
  tp53_upper <- tp53_quantiles[[4]]
  tp53_lower <- tp53_quantiles[[2]]
  telomerase_upper <- telomerase_quantiles[[4]]
  telomerase_lower <- telomerase_quantiles[[2]]
  
  meta_indep <- meta_indep %>%
    dplyr::mutate(
      # use 0.5 as a cutoff for tp53 score and telomerase score
      tp53_0.5 = dplyr::case_when(
        tp53_score > 0.5  ~ "tp53_high",
        tp53_score <= 0.5 ~ "tp53_low", 
        TRUE ~ NA_character_),
      tel_0.5 = dplyr::case_when(
        telomerase_score > 0.5  ~ "telomerase_high",
        telomerase_score <= 0.5 ~ "telomerase_low",
        TRUE ~ NA_character_),
      
      # for strict ones, 25th and 75th quantile is used
      tp53_strict = dplyr::case_when(
        tp53_score >= tp53_upper ~ "tp53_high",
        tp53_score <= tp53_lower ~ "tp53_low",
        tp53_score > tp53_lower & tp53_score < tp53_upper ~ "tp53_mid",
        TRUE ~ NA_character_),
      tel_strict = dplyr::case_when(
        telomerase_score > telomerase_upper ~ "telomerase_high",
        telomerase_score <= telomerase_lower ~ "telomerase_low",
        telomerase_score > telomerase_lower & telomerase_score <= telomerase_upper ~ "telomerase_mid",
        TRUE ~ NA_character_)
    )
      
  
  meta_indep <- meta_indep %>% 
    # for pheno_0.5, this is a column where combined info from tp53_0.5 and tel_0.5 
    dplyr::mutate(
      pheno_0.5 = ifelse(is.na(tp53_0.5) | is.na(tel_0.5), NA, paste(tp53_0.5, tel_0.5))
    ) 
    
  meta_indep <- meta_indep %>% 
    # for pheno_strict, this is a column combining info about whether the sample falls into any extreme category 
    # or locates in the middle for either tp53 or telomerase - basically tp53_strict and tel_strict combined
    dplyr::mutate(
      pheno_strict = ifelse(is.na(tp53_strict) | is.na(tel_strict), NA, paste(tp53_strict, tel_strict))
    )
  
  meta_indep <- meta_indep %>% 
    # this is a column only retains info from pheno_strict when both tp53 and telomerase do not fall into mid category
    dplyr::mutate(
      pheno_extremes = dplyr::case_when(
        !grepl("mid", pheno_strict) ~ pheno_strict, 
        TRUE ~ NA_character_)
    )
  
  ### get the results out 
  return(meta_indep)
}
