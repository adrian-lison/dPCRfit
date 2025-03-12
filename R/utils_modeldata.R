#' Create a model_component object
#'
#' @param name A character string representing the name of the model component.
#' @param .e An expression that defines the model component through changes to
#'   the modeldata object.
#' @param args Arguments to be passed to the model component. By default, this
#'   is taken from the calling function.
#'
#' @return An S3 object of class "model_component".
#' @export
#' @keywords internal
model_component <- function(name, .e, args = as.list(parent.frame())) {
  if (!is.character(name) || length(name) != 1) {
    cli::cli_abort("`name` must be a single character string")
  }
  if (!is.list(args)) {
    cli::cli_abort("`args` must be a list")
  }

  func <- function() {}
  formals(func) <- c(args, list(modeldata = NULL))
  body(func) <- substitute(.e)

  component <- do.call(purrr::partial, c(list(func), args))

  structure(
    list(name = name, args = args, component = component),
    class = "model_component"
  )
}

#' @export
setClass("model_component")

#' Print method for model_component
#'
#' @param x A model_component object.
#' @param ... Additional arguments passed to print.
#' @export
#' @keywords internal
print.model_component <- function(x, ...) {
  cat("<dPCRfit model_component>", "\n")
  cat("Name:", x$name, "\n")
  show_args <- sapply(x$args, function(a) if(is.data.frame(a)) "<data.frame>" else a)
  cat("Arguments:\n", paste(names(x$args), show_args, sep = " = " , collapse = "\n "), "\n")
}

#' Construct an unspecified model
#'
#' @return A `modeldata` object containing data and specifications of the model
#'   to be fitted. Can be combined with other model components to
#'   add data and model specifications.
#'
#' @return The `modeldata` object also includes information about parameter
#'   initialization (`.init`), meta data (`.metainfo`), and checks to be
#'   performed before model fitting (`.checks`).
#' @export
modeldata_init <- function() {
  modeldata <- list()
  modeldata$.inputs <- list()
  modeldata$.init <- list()
  modeldata$.metainfo <- list()
  modeldata$.checks <- list()

  class(modeldata) <- "modeldata"
  return(modeldata)
}

#' @export
setClass("modeldata")

#' @export
`+.modeldata` <- function(modeldata, x) {
  if (!inherits(x, "model_component")) {
    cli::cli_abort("You can only add model_component objects to the modeldata.")
  }
  return(x$component(modeldata = modeldata))
}
