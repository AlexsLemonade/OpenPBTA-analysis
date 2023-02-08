`%>%` <- dplyr::`%>%`


#' Filter an expression data frame or matrix by median cancer group tumor purity
#'
#' @param expression_data A data frame or matrix containing expression data. Columns are expected to be genes, 
#'   and rows are expected to be samples. 
#' @param bs_ids_file Path to TSV file which contains, at least, a `Kids_First_Biospecimen_ID`
#'   column that contains only IDs that pass a given threshold. 
#'
#' @return An updated `expression_data` object that contains expression data only for
#'  samples that pass the threshold
filter_expression_by_tumor_purity <- function(expression_data, 
                                              bs_ids_file){
  
  # Define vector of IDs to retain
  bs_ids <- readr::read_tsv(bs_ids_file) %>%
    dplyr::pull(Kids_First_Biospecimen_ID) %>%
    unique()
  
  # Subset the provided variable to these IDs
  # To handle that the input can either be a data frame or matrix, we first
  #  subset `bs_ids` to only the intersection with the input's names.
  #  This allows us to use the same code to filter on either input object structure.
  bs_ids <- intersect(bs_ids, colnames(expression_data))
  
  # Warn if some ids are missing from the expression data.
  missing_bs_ids <- setdiff(bs_ids, colnames(expression_data))
  if ( length(missing_bs_ids) > 0 ) {
    collapsed <- paste(missing_bs_ids, collapse = "\n")
      warning(
      glue::glue("The following ids are missing from the expression data: 
        {collapsed}")
      )
  }
  
  # Subset the `expression_data` variable to those IDs
  expression_data <- expression_data[, bs_ids]
  
  # Return the updated `expression_data` object
  return(expression_data)
}
