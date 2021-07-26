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
# - Only functions that are not nested in other expressions can be imported,
#   e.g. functions nested in other functions cannot be imported
# - Only functions defined by `<-` or `assign` can be imported, e.g. unnamed
#   functioins cannot be imported
# - The import_function raises an error if the imported function is not found in
#   the file, or the imported function is defined multiple times in the
#   source_code_R_file
# - If the imported function relies on external variables, the variables will
#   not be immediately available after this import_function returns. Users could
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
  # envir = parent.frame() makes the environment(imported_function) to be the
  # same as the import function being called. Even though the eval documentation
  # says the default envir parameter is parent.frame(), leaving envir =
  # parent.frame() in the eval call will surprisingly make the
  # environment(imported_function) to be the same as the environment of the
  # import_function call. Maybe this is caused by the place where parent.frame()
  # is evaluated? If specified, parent.frame() is evaluated in the calling
  # function; if not specified, parent.frame() is evaluated in the eval call
  # environment?

  # Examples executed in terminal R, in order to avoid RStudio customizations.
  # > print(environment())
  # <environment: R_GlobalEnv>
  # >
  # > foo <- function() {
  # +   print(parent.frame())
  # +   print(environment())
  # +   return(eval(quote(function() { return(2) })))
  # + }
  # >
  # > bar <- foo()
  # <environment: R_GlobalEnv>
  # <environment: 0x5645ead67670>
  # > print(environment(bar))
  # <environment: 0x5645ead67670>
  # >
  # > baz <- function() {
  # +   print(parent.frame())
  # +   print(environment())
  # +   return(eval(quote(function() { return(2) }),
  # +               envir = parent.frame()))
  # + }
  # >
  # > qux <- baz()
  # <environment: R_GlobalEnv>
  # <environment: 0x5645ead5fec0>
  # > print(environment(qux))
  # <environment: R_GlobalEnv>

  # Only evaluate the function() {} part
  return(eval(function_def_exprs[[1]][[3]], envir = parent.frame()))
}
