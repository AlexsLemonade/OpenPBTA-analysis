############################################################################################
#    Functions for performing ANOVA and Tukey HSD tests on hallmark pathway GSVA scores    #
#                                                                                          #
#    Stephanie J. Spielman, 2020                                                           #
############################################################################################

library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(broom)


# Magrittr pipe
`%>%` <- dplyr::`%>%`

### Check tidyr version for nesting code to ensure it used as <1.0.0 ####
tidyr_version <- package_version( packageVersion("tidyr") )
if (tidyr_version$major >= 1){
  nest   <- nest_legacy
  unnest <- unnest_legacy
}




perform_anova <- function(df, predictor_variable_string)
{
  ### Defines a function for performing ANOVA to detect signficance GSVA scores according to a given predictor variable
  ### Specifically to be used within `compare_pathways()`
  formula <- as.formula(paste("gsea_score ~", predictor_variable_string))
  aov(formula, data = df)
}

gsva_anova_tukey <- function(df, predictor_variable, library_type, significance_threshold)
{
  ### Defines a function to compare GSVA scores across a given predictor variable, _for each pathway (hallmark)_.
  ### Returns a named list with two items:
  ##### 1. "anova" is a tibble of all ANOVAs performed, for each hallmark
  ##### 2. "tukey" is a tibble of all Tukey HSD tests performed, comparing each level of `predictor_variable` for each hallmark
  
  ### Input arguments:
  ### 1. `df`:   the data frame containing scores (`gsea_score`) and assumed to contain the specified `predictor_variable` and a column `data_type` referring to either "stranded" or "polyA" libraries
  ### 2. `predictor_variable`: A _string_ indicating the predictor variable, such as `short_histology`, `disease_type_new`, etc.
  ### 3. `library_type`:  Which RNAseq library to analyze? Specify either `polyA` or `stranded`
  ### 4. `significance_threshold`:  What threshold is used for significance (aka alpha, FPR)?

  predictor_variable        <- enquo(predictor_variable)  
  predictor_variable_string <- quo_name(predictor_variable)
  library_type <- str_to_lower(library_type)
  ######## Sanity check #########
  if (!(predictor_variable_string %in% colnames(df))) stop("ERROR: The provided `predictor_variable` is not found in the data frame.")
  if (!("hallmark_name" %in% colnames(df))) stop("ERROR: The needed column `hallmark_name` is not found in the data frame.")
  if (!("gsea_score" %in% colnames(df))) stop("ERROR: The needed column `gsea_score` is not found in the data frame.")
  if (!(library_type %in% unique(df$data_type))) stop(paste("ERROR: The provided `library_type` must be one of:", 
                                                         paste(unique(df$data_type), collapse= " or "),
                                                         "."))

  ### Determine how many ANOVAs we will be doing
  df %>% 
    dplyr::select(hallmark_name, data_type) %>%
    distinct() %>% 
    count(data_type) %>%
    filter(data_type == library_type) %>%
    pull(n) -> number_of_tests 
  
  print(number_of_tests)
  ############### Perform modeling, including ANOVAs and Tukey across hallmarks
  df %>%
    filter(data_type == library_type) %>%
    group_by(hallmark_name) %>%
    nest() %>%
    mutate(anova_fit = pmap(list(data, predictor_variable_string), perform_anova),
           tukey_fit = map(anova_fit, TukeyHSD)) -> fitted_results
  
  ############## Tidy the ANOVA fits
  fitted_results %>%
    dplyr::select(hallmark_name, anova_fit) %>%
    mutate(anova_tidy = map(anova_fit, broom::tidy)) %>%
    dplyr::select(-anova_fit) %>%
    unnest() %>%
    ungroup() %>%
    filter(term != "Residuals") %>%
    dplyr::select(-df, -sumsq, -meansq, -term) %>%
    rename(anova_f_statistic = statistic,
           anova_p_value     = p.value) %>%
    mutate(anova_p_value_bonferroni = anova_p_value * number_of_tests,
           anova_p_value_bonferroni = ifelse(anova_p_value_bonferroni >= 1, 1, anova_p_value_bonferroni),
           significant_anova = ifelse(anova_p_value_bonferroni <= significance_threshold, TRUE, FALSE)) -> final_anova_results
    
  ############## Tidy the Tukey tests
  fitted_results %>%
    dplyr::select(hallmark_name, tukey_fit) %>%
    mutate(tukey_tidy = map(tukey_fit, broom::tidy)) %>%
    dplyr::select(-tukey_fit) %>%
    unnest() %>%
    dplyr::select(hallmark_name, comparison, estimate, adj.p.value) %>%
    dplyr::rename(pathway_score_difference = estimate,
                  tukey_p_value            = adj.p.value) %>%
    mutate(bonferroni_pvalue = tukey_p_value * number_of_tests,
           bonferroni_pvalue = ifelse(bonferroni_pvalue >= 1, 1, bonferroni_pvalue), 
           significant_tukey = tukey_p_value <= significance_threshold,
           significant_tukey_bonf = bonferroni_pvalue <= significance_threshold) -> final_tukey_results

  
  return(list("anova" = final_anova_results, "tukey" = final_tukey_results))
}
