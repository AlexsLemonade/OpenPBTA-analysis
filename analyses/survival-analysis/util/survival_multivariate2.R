# Function to fit, visualize, and export a cox regression. 
# A table of results is exported, but plots are not exported
# This is specifically useful for **multivariate regressions**

fit_plot_save_cox_multivariate2 <- function(df,    # data frame with the data (assumed NOT recoded)
                                           terms, # string of covariates
                                           covariate_dict, # dictionary of covariates to readable names
                                           ref_dict, # dictionary of covariate reference levels
                                           output_dir, 
                                           filename = NULL,
                                           plot_dir = NULL) {
  # df: A data frame that contains columns:
  #   - `OS_status` with character values "0" and "1"
  #   - `OS_years` giving survival time in years
  # terms: A string providing RHS of model equation, using terms that match column names in `df`
  # covariate_dict: dictionary mapping column names to readable names for table/plot
    # Example: covariate_dict = c(age="Age at Dx", sex="Sex")
  # ref_dict: dictionary mapping reference level for any/all covariates
    # Example: ref_dict = c(hgg_group="non-HGAT")
  # output_dir: Output directory to save model result table in. 
  # filename: (Optional) Output file name for model result table. If nothing is given, it will be derived from model terms.
  # plot_dir: (Optional) Output folder for forest plot.
  # custom_x_breaks: (Optional) Custom x-axis breaks for HR plot (Default: c(0.25, 0.5, 0.75, 1, 1.5, 2)) - add this later
  
  # Define and fit the model
  fit_result <- survivalAnalysis::analyse_multivariate(df,
                                                vars(time = OS_years, status = OS_status),
                                                covariates = terms,
                                                covariate_name_dict = covariate_dict,
                                                reference_level_dict = ref_dict)
  

  # Save results table fit
  if (is.null(filename)) {
    # do not include * in file names
    filename <- paste0("SurvivalAnalysis_cox_reg_results_per_", 
                       stringr::str_replace_all(cov, "\\*", "-by-"),
                       ".tsv")
  } 
  
  readr::write_tsv(fit_result$summaryAsFrame,
                   file.path(output_dir, filename)
  ) 
  
  # Create visualization
  if(is.null(plot_dir)) {
  }
  forest_coxph2 <- survivalAnalysis::forest_plot(fit_result,
                                                 factor_labeller = covariate_dict,
                                                 endpoint_labeller = c(OS_years="OS"),
                                                 orderer = ~order(HR),
                                                 labels_displayed = c("endpoint", "factor", "n"),
                                                 ggtheme = ggplot2::theme_bw(base_size = 8),
                                                 relative_widths = c(1.5, 2, 1))
                                                 #HR_x_breaks <- c(0.25, 0.5, 0.75, 1, 1.5, 2))
  
  # Print the tidied fit and the plot
  print(fit_result$summaryAsFrame)
  print(forest_coxph2)
  
  save_pdf(plot = forest_coxph2, 
           folder = plot_dir, 
           fileBaseName = paste0("HR_", cg, "_", cov))
  
}
