# Functions for conducting survival analyses
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################

survival_analysis <- function(metadata, ind_var, test = "kap.meier", data = NULL) {
  # Given the overall survival information, and an independent variable will 
  # run survival analysis and return a list with 1) the model fit object and 
  # 2) a summary table. 
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
  #   data: If the data for the independent variable needed for the test is not
  #         in metadata, add it here. This assumes it will be a data.frame with 
  #         either a "Tumor_Sample_Barcode" or "Kids_First_Biospecimen_ID" column 
  #         which to inner_join with metadata by. 
  #
  # Returns:
  # A list with two objects: the model fit object and the summary table
  
  # List the tests
  supported_tests <- c("kap.meier", "log.rank", "cox.reg")
  
  # Check that it is a supported test
  if (!(test %in% supported_tests)){
    stop(paste0(test, "is not a supported test."))
  }

  # Format the OS_status variable for the survival model functions
  metadata <- readr::read_tsv(metadata_file)
  metadata$OS_status <- factor(metadata$OS_status, levels = c("LIVING", "DECEASED"))
  metadata$OS_status_num <- as.numeric(metadata$OS_status)
  
  # If other data has been supplied, attempt to join it
  if (!is.null(data)) {
    # Make it so the column names are automatically matching no matter what
    colnames(data) <- gsub("Tumor_Sample_Barcode", "Kids_First_Biospecimen_ID", colnames(data))
    
    # Join this data to the metadata
    metadata <- dplyr::inner_join(metadata, data, 
                                  by = "Kids_First_Biospecimen_ID")
  }

  # Extract independent variables
  ind_var_df <- metadata %>% 
    dplyr::select(colnames(metadata) %in% eval(parse(text = ind_var)))
  
  # For the model need a plus sign for separating multiple independent variables
  if(length(ind_var) > 1) {
    ind_var <- paste0(ind_var, collapse = "+")
  }
  
  # Piece together a model
  model <- paste("survival::Surv(metadata$OS_days, metadata$OS_status_num)", 
                 ind_var, 
                 sep = " ~ ")
  
  # Print out what the model is
  message(paste("Testing model:", model))
  
  # Run the appropriate test
  if (test == "kap.meier") {
    # Make the model
    fit <- survival::survfit(
      survival::Surv(metadata$OS_days, metadata$OS_status_num) ~ eval(parse(ind_var)),
      data = metadata 
    )
  }
  if (test == "log.rank") {
    # Make the model
    fit <- survival::survdiff(
      survival::Surv(metadata$OS_days, metadata$OS_status_num) ~ eval(parse(ind_var)),
      data = metadata
    )
  }
  if (test == "cox.reg") {
    # Make the model
    fit <- survival::coxph(
      survival::Surv(metadata$OS_days, metadata$OS_status_num) ~ eval(parse(ind_var)),
      data = metadata
    )
  }
  
  # Tidy up the model object with broom
  table <- broom::tidy(fit) 
  
  # Return both the fit object and the table
  return(list(fit, table))
}
