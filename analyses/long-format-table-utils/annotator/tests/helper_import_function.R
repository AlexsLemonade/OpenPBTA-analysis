# Import a single function from an R file without running the whole file
#
# Args:
# - source_code_R_file: path to the source code R file that contains the
#   definition of the function to be imported
# - function_name: the name of the function to be imported
#
# Returns the imported function
#
# NOTES:
# - Only functions defined by `<-`` or `assign` can be imported
# - The function returns an error if the function is not found in the file, or
#   the function is defined multiple times in the file
# - If the imported function relies on external variables, the variables will
#   not be immediately available after this function completes. Users could
#   import dependency functions and define dependency variables at a later
#   point.
#
# Adapted from @MrFlick's answer at https://stackoverflow.com/a/29132294/4638182
import_function <- function(source_code_R_file, function_name) {
  src_exprs <- parse(source_code_R_file)
  function_def_exprs <- purrr::keep(
    src_exprs,
    function(x) {
      # x[[1]] has class name, which is not identical to "<-"
      if(x[[1]] == "<-" || x[[1]] == "assign") {
        # can only have [[2]] and [[3]] if x[[1]] == "<-"
        if(x[[2]] == function_name && x[[3]][[1]] == "function") {
          return(TRUE)
        }
      }
      return(FALSE)
    }
  )
  if (identical(length(function_def_exprs), as.integer(0))) {
    stop(paste0(function_name, " is not found in ", source_code_R_file, "."))
  } else if (!identical(length(function_def_exprs), as.integer(1))) {
    stop(paste0(function_name, " is found multiple times in ",
                source_code_R_file, "."))
  }
  # only evaluate the function() {} part
  return(eval(function_def_exprs[[1]][[3]]))
}
