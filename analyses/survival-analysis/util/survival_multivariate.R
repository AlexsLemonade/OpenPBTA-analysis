# Function to fit, visualize, and export a cox regression. 
# A table of results is exported, but plots are not exported
# This is specifically useful for **multivariate regressions**
fit_plot_save_cox_multivariate <- function(df,    # data frame with the data (assumed NOT recoded)
                                           terms, # string giving predictors
                                           output_dir, 
                                           filename = NULL, 
                                           plot_title = NULL) {
  # df: A data frame that contains columns:
  #   - `OS_status` with character values "LIVING" and "DECEASED"
  #   - `OS_years` giving survival time in years
  # terms: A string providing RHS of model equation, using terms that match column names in `df`
  # output_dir: Output directory to save model result table in. 
  # filename: (Optional) Output file name for model result table. If nothing is given, it will be derived from model terms.
  # plot_title: (Optional) A title for the plot that will be printed. If nothing is given, it will be derived from model terms.
  
  # Recode OS_Status
  df <- df %>%
    mutate(OS_status = ifelse(OS_status == "LIVING", 0, 1))
  
  # Define and fit the model
  model <- paste0("survival::Surv(time = OS_years, event = OS_status) ~ ", terms)
  fit <- survival::coxph(
    formula(model),
    data = df
  )
  
  tidy_fit <- broom::tidy(fit)
  
  # Save tidied fit
  if (is.null(filename)) {
    # do not include * in file names
    filename <- paste0("cox_reg_results_per_", 
                       stringr::str_replace_all(terms, "\\*", "-by-"),
                       ".tsv")
  } 
  
  readr::write_tsv(tidy_fit,
                   file.path(results_dir, filename)
  ) 
  
  # Create visualization
  if(is.null(plot_title)) {
    plot_title <- terms
  }
  forest_coxph <- survminer::ggforest(fit, data = df, main = paste0("Hazard ratio for ", plot_title))
  
  # Print the tidied fit and the plot
  print(tidy_fit)
  print(forest_coxph)

  
}
