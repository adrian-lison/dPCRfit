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

suppress_warnings <- function(.expr, .f, ...) {
  eval.parent(substitute(
    withCallingHandlers(.expr, warning = function(w) {
      cm <- conditionMessage(w)
      cond <- any(sapply(.f_list, function(.f) {
        if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
      }
      ))
      if (cond) {
        invokeRestart("muffleWarning")
      }
    })
  ))
}

suppress_messages <- function(.expr, .f_list, ...) {
  eval.parent(substitute(
    withCallingHandlers(
      .expr,
      message = function(w) {
        cm <- conditionMessage(w)
        cond <- any(sapply(.f_list, function(.f) {
          if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
        }
        ))
        if (cond) {
          invokeRestart("muffleMessage")
        }
      }
    )
  ))
}

suppress_messages_warnings <- function(.expr, .f_list, ...) {
  eval.parent(substitute(
    withCallingHandlers(
      .expr,
      message = function(w) {
        cm <- conditionMessage(w)
        cond <- any(sapply(.f_list, function(.f) {
          if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
        }
        ))
        if (cond) {
          invokeRestart("muffleMessage")
        }
      },
      warning = function(w) {
        cm <- conditionMessage(w)
        cond <- any(sapply(.f_list, function(.f) {
          if (is.character(.f)) grepl(.f, cm) else rlang::as_function(.f)(cm, ...)
        }
        ))
        if (cond) {
          invokeRestart("muffleWarning")
        }
      }
    )
  ))
}

#' Semantic sugar to wrap an abort statement in a function
#'
#' This can be used to provide abort statements as a direct argument to
#' tryCatch without wrapping them in a function call
#' @keywords internal
abort_f <- function(...) {
  return(function(err) {
    cli::cli_abort(
      call = NULL, .trace_bottom = sys.frames()[[sys.nframe() - 4]], ...
    )
  })
}

#' Collect warnings arising during computation of a result
#'
#' @param expr Expression with the computation to run
#'
#' @return A `list`, with two elements: value and warnings
#' @keywords internal
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}

#' @keywords internal
format_table <- function(df) {
  # Ensure df is a data frame
  if (!is.data.frame(df)) stop("Input must be a data frame")

  # Convert all columns to character for consistent formatting (numeric columns with max 2 decimals)
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      format(col, digits = 4)
    } else {
      as.character(col)
    }
  })

  # Determine max width per column
  col_widths <- sapply(df, function(col) max(nchar(col), na.rm = TRUE))
  col_widths <- pmax(col_widths, nchar(names(df)))  # Ensure headers are included

  # Function to format each row with proper spacing
  format_row <- function(row) {
    paste(mapply(sprintf, paste0("%-", col_widths, "s"), row), collapse = "  ")
  }

  # Format header and rows
  header <- format_row(names(df))
  rows <- apply(df, 1, format_row)

  # Combine into a single text output
  table_text <- paste(c(header, rows), collapse = "\n")

  return(table_text)
}
