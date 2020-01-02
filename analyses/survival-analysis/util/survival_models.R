# Functions for conducting survival analyses
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Install the packages we need
if (!("survminer" %in% installed.packages())) {
  install.packages("survminer")
}

# Attach this package
library(survminer)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

survival_analysis <- function(metadata,
                              ind_var,
                              test = "kap.meier",
                              ind_data = NULL,
                              metadata_sample_col = "Kids_First_Biospecimen_ID",
                              ind_data_sample_col = NULL,
                              os_days_col = "OS_days",
                              os_status_col = "OS_status") {
  # Given the overall survival information, and an independent variable will
  # run survival analysis and return a list with 1) the model fit object and
  # 2) a summary table 3) the data used.
  #
  # Args:
  #   metadata: a data.frame that contains columns OS_status and OS_days to use
  #             for the survival model. This also assumes "LIVING" and "DECEASED"
  #             are the two statuses. THis will be converted to a numeric variable
  #             for use with `survival` R functions. Samples with NAs are dropped.
  #   ind_var: a character string noting the name of the independent variable to
  #            test as a predictor for survival.
  #   test: a character string noting which test to use. Supported choices:
  #         "kap.meier", "log.rank", "cox.reg"
  #   ind_data: If the data for the independent variable needed for the test is not
  #         in metadata, add it here. This assumes it will be a data.frame with
  #         either a "Tumor_Sample_Barcode", "Kids_First_Biospecimen_ID", or
  #         "sample" column which to inner_join with metadata by.
  #   ind_data_sample_col : A character string that states the name of the column in
  #                `ind_data` that contains the sample IDs that match to the
  #                metadata.
  #   metadata_sample_col : A character string that states the name of the column in
  #                `metadata` that contains the sample IDs that match to the
  #                ind_data. Default is "Kids_First_Biospecimen_ID".
  #   os_days_col : A character string that states the name of the column in
  #                `metadata` that contains the overall survival days information.
  #                 Default is "OS_days".
  #   os_status_col : A character string that states the name of the column in
  #                `metadata` that contains the overall survival status information
  #                 The data in this column will be converted to numeric if it
  #                 is not already.  Default is "OS_status".
  #
  # Returns:
  # A list with three objects: 1) the original model fit object,
  #                            2) the summary table
  #                            3) the original data.frame

  ####################### Check the options supplied ###########################
  # List the tests
  supported_tests <- c("kap.meier", "log.rank", "cox.reg")

  # Check that it is a supported test
  if (!(test %in% supported_tests)) {
    stop(paste(
      test, "is not a supported test. Please specify one of the following:",
      paste(supported_tests, collapse = ", ")
    ))
  }

  # Don't continue if there's no variable specified
  if (is.null(ind_var)) {
    stop("No variable has been supplied to test with using the `ind_var` argument. 
         Stopping.")
  }

  # List the columns we need
  needed_cols <- c(os_status_col, os_days_col)

  # Get logical vector indicating which are in metadata
  found_cols <- (needed_cols %in% colnames(metadata))

  # If not all the columns are found, stop
  if (!all(found_cols)) {
    stop(cat(
      "The following column names specified for overall survival information: \n",
      paste(needed_cols[which(found_cols)], collapse = "\n"),
      "\n ...were not found in the specified metadata data.frame.",
      "Check your `os_status_col` and `os_days_col` arguments."
    ))
  }

  ######################## Set up the survival variables #######################
  # Pull out this data
  os_status <- metadata %>%
    dplyr::pull({{os_status_col}})

  # Reformat as a numeric variable where 1 = LIVING and 2 = DECEASED.
  metadata <- metadata %>%
    dplyr::mutate(
      !!os_status_col := as.numeric(
        factor(os_status, levels = c("LIVING", "DECEASED")))
      )

  ############################ Set up the ind data #############################
  # If other ind_data has been supplied, attempt to join it to metadata
  if (!is.null(ind_data)) {

    # Check that a sample column was specified
    if (is.null(ind_data_sample_col)) {
      stop("A `ind_data` data.frame was specified but the column with the sample information 
           was not specified using `ind_data_sample_col`.")
    }

    # List the columns we need
    found_cols <- c(
      ind_data_sample_col %in% colnames(ind_data),
      ind_var %in% colnames(ind_data),
      metadata_sample_col %in% colnames(metadata)
    )

    # If not all the columns are found, stop
    if (!all(found_cols)) {
      stop(cat(
        "The following column names specified: \n",
        paste(needed_cols[which(found_cols)], collapse = "\n"),
        "\n ...were not found in the specified data.frames.",
        "Check your `metadata_sample_col`, `ind_data_sample_col`, and `ind_var` arguments."
      ))
    }

    # Set up the enquosures
    ind_data_sample_col <- enquo(ind_data_sample_col)
    metadata_sample_col <- enquo(metadata_sample_col)
    
    # Set up the sample by thing
    sample_by <- set_names(quo_name(ind_data_sample_col), quo_name(metadata_sample_col))

    # Join this ind_data to the metadata
    metadata <- dplyr::inner_join(metadata, ind_data,
      by = sample_by
    )
  } else {
    # Set up
    metadata_sample_col <- enquo(metadata_sample_col)
    
    # Look for the independent variable columns in the metadata
    found_cols <- (ind_var %in% colnames(metadata))

    # If not all the columns are found, stop
    if (!all(found_cols)) {
      stop(cat(
        "The following column names specified for the independent variable(s): \n",
        ind_var[which(found_cols)],
        "Check your `ind_var` argument."
      ))
    }
  }

  ############################ Set up model and data ###########################
  # Extract independent variables
  ind_var_df <- metadata %>%
    dplyr::select(!!metadata_sample_col, {{os_days_col}}, !!os_status_col, ind_var)

  # For the model need a plus sign for separating multiple independent variables
  ind_var <- paste0(ind_var, collapse = "+")

  # Piece together a model
  model <- paste0(
    "survival::Surv(",
    os_days_col, ", ", os_status_col, ") ~ ",
    ind_var
  )

  ################################# Do the test! ###############################
  # Print out what the model is
  message(paste("Testing model:", model, "with", test))

  # Run the appropriate test
  if (test == "kap.meier") {
    # Make the model
    fit <- survival::survfit(
      formula(model),
      data = ind_var_df
    )
  }
  if (test == "log.rank") {
    # Make the model
    fit <- survival::survdiff(
      formula(model),
      data = ind_var_df
    )
    # Obtain p value for Chi-Squared stat
    fit$p.value <- pchisq(fit$chisq, df = 1, lower = FALSE)
  }
  if (test == "cox.reg") {
    # Make the model
    fit <- survival::coxph(
      formula(model),
      data = ind_var_df
    )
  }
  # Tidy up the model object with broom
  table <- broom::tidy(fit)

  # Restore the model in this slot so the ggsurvplot function can find it
  fit$call$formula <- formula(model)

  # Return both the fit object and the table
  return(list(model = fit, table = table, original_data = ind_var_df))
}
