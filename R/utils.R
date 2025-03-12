#' Print link to help page of function in CLI
#'
#' @param function_name Name of the function
#' @param label Optional label for the link, if NULL (default), the function
#'   name is used as label
#'
#' @return A string with a link to the help page of the function
#' @keywords internal
cli_help <- function(function_name, label = NULL) {
  function_name <- paste0(function_name, "()")
  if (is.null(label)) {
    label <- function_name
  }
  return(paste0(
    "{.help [",
    label,
    "](dPCRfit::",
    function_name,
    ")}"
  ))
}
