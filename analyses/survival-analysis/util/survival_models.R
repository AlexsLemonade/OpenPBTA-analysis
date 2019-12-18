# Functions for conducting survival analyses
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
# Install the packages we need
if (!("survival" %in% installed.packages())) {
  install.packages("survival")
}

if (!("survminer" %in% installed.packages())) {
  install.packages("survminer")
}

# Attach this package
library(survminer)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

survival_analysis <- function(metadata, ind_var, test = "kap.meier", ind_data = NULL) {
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
  #
  # Returns:
  # A list with three objects: 1) the original model fit object,
  #                            2) the summary table
  #                            3) the original data.frame

  # List the tests
  supported_tests <- c("kap.meier", "log.rank", "cox.reg")

  # Check that it is a supported test
  if (!(test %in% supported_tests)) {
    stop(paste0(test, "is not a supported test."))
  }

  # Format the OS_status variable for the survival model functions
  metadata$OS_status <- factor(metadata$OS_status, levels = c("LIVING", "DECEASED"))
  metadata$OS_status_num <- as.numeric(metadata$OS_status)

  # If other ind_data has been supplied, attempt to join it
  if (!is.null(ind_data)) {
    # Make it so the column names are automatically matching no matter what
    colnames(ind_data) <- gsub("Tumor_Sample_Barcode|sample", "Kids_First_Biospecimen_ID", colnames(ind_data))

    # Join this ind_data to the metadata
    metadata <- dplyr::inner_join(metadata, ind_data,
      by = "Kids_First_Biospecimen_ID"
    )
  }

  # Extract independent variables
  ind_var_df <- metadata %>%
    dplyr::select("Kids_First_Biospecimen_ID", "OS_status_num", "OS_days", ind_var)

  # For the model need a plus sign for separating multiple independent variables
  ind_var <- paste0(ind_var, collapse = "+")

  # Piece together a model
  model <- paste("survival::Surv(OS_days, OS_status_num)",
    ind_var,
    sep = " ~ "
  )
  
  # We'll need to grab this later for ggsurvplot
  assign("model", model, envir = .GlobalEnv)
  
  # Print out what the model is
  message(paste("Testing model:", model))

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

  # Return both the fit object and the table
  return(list(model = fit, table = table, original_data = ind_var_df))
}
