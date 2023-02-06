


#' Filter an expression data frame or matrix by median cancer group tumor purity
#'
#' @param expression_data A data frame or matrix containing expression data. Columns are expected to be genes, 
#'   and rows are expected to be samples. 
#' @param tumor_purity_file Path to tumor purity metadata file TSV which contains, at least, a `Kids_First_Biospecimen_ID`
#'   column that has been filtered to IDs that pass a give threhsold. The default file contains biospecimens that pass 
#'   the cancer-group-specific median threshold.
#'
#' @return An updated `expression_data` object that contains expression data only for
#'  samples that pass the threshold
filter_expression_by_tumor_purity <- function(expression_data, 
                                              tumor_purity_file = 
                                                file.path(
                                                  rprojroot::find_root(rprojroot::has_dir(".git")),
                                                  "analyses",
                                                  "tumor-purity-exploration",
                                                  "results",
                                                  "thresholded_rna_stranded_same-extraction.tsv"
                                                )
                                               ){
  
  # Define vector of IDs to retain
  bs_ids <- readr::read_tsv(tumor_purity_file) %>%
    dplyr::pull(Kids_First_Biospecimen_ID) %>%
    unique()
  
  # Subset the provided variable to these IDs
  # To handle that the input can either be a data frame or matrix, we first
  #  subset `bs_ids` to only the intersection with the input's names.
  #  This allows us to use the same code to filter on either input object structure.
  bs_ids <- intersect(bs_ids, colnames(expression_data))
  
  missing_bs_ids <- setdiff(bs_ids, colnames(expression_data))
  if ( length(missing_bs_ids) > 0 ) {
     warn(glue::glue("{missing_bs_ids} is missing from expression data"))
  } 
  
  # Subset the `expression_data` variable to those IDs
  expression_data <- expression_data[, bs_ids]
  
  # Return the updated `expression_data` object
  return(expression_data)
}
