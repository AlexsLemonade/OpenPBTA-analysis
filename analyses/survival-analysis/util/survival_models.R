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
  # Run survival analysis. 
  #
  # Args:
  # metadata:
  # ind_var:
  # test:
  # data:
  #
  # Returns:
  # A plot or table depending on the test specified
  
  
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
  
  return(table)
}
