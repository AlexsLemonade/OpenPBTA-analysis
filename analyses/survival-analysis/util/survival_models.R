# Functions for conducting survival analyses
#
# C. Savonen for ALSF - CCDL
#
# 2019
#
################################################################################
########################### Setting Up Functions ###############################
################################################################################

survival_analysis <- function(metadata, ind_var, data = NULL) {
  # Run survival analysis. 
  #
  # Args:
  # metadata:
  # ind_var:
  # data:
  #
  # Returns:
  # A plot or table
  #
  # Format the OS_status variable for the survival model functions
  metadata <- readr::read_tsv(metadata_file)
  metadata$OS_status <- factor(metadata$OS_status, levels = c("LIVING", "DECEASED"))
  metadata$OS_status_num <- as.numeric(metadata$OS_status)
  
  # If other data has been supplied, attempt to join it
  if (!is.null(data)) {
    metadata <- dplyr::inner_join(metadata, data, 
                                  by = c("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode"))
  }

  # Extract independent variables
  ind_var_df <- metadata %>% 
    dplyr::select(colnames(metadata) %in% eval(ind_var))
  
  model <- paste("survival::Surv(metadata$OS_days, metadata$OS_status_num)", 
                 paste(ind_var, collapse = " + "), 
                 sep = " ~ ")
  
  message(paste("Testing model:", model))
  
  # The new line of code
  model <- eval(model)
  
  # Make the model
  kap_fit <- survival::survfit(
    eval(model),
    data = germline_metadata 
  )
  
  # Print out a cleaned version of this data to look at it
  broom::tidy(kap_fit)
  
}
