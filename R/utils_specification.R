select_required_cols <- function(df, required_cols, col_names){

  col_names <- col_names[!sapply(required_cols,is.null)]
  required_cols <- purrr::list_c(required_cols)

  if (!all(required_cols %in% names(df))) {
    cli::cli_abort(
      c(paste(
        "The following columns must be present",
        "in the provided measurements `data.frame`:",
        paste(required_cols, collapse = ", ")
      ),
      paste("Please adjust the `data.frame` or specify the right column",
            "names via the `_col` arguments of this function.")
    )
    )
  }

  df = as.data.table(df)[, .SD, .SDcols = required_cols]
  df <- setnames(df, old = required_cols, new = col_names)
  return(df)
}

set_distribution <- function(distribution) {
  available_distributions <- c(
    "gamma" = 0,
    "log-normal" = 1,
    "truncated normal" = 2,
    "normal" = 3
  )

  distribution_aliases <- c(
    "lnorm" = "log-normal",
    "lognormal" = "log-normal",
    "truncated_normal" = "truncated normal",
    "truncated-normal" = "truncated normal",
    "norm" = "normal"
  )

  distribution <- stringr::str_to_lower(distribution)
  distribution <- do.call(
    switch, c(distribution, as.list(c(distribution_aliases,distribution)))
  )

  if (!distribution %in% names(available_distributions)) {
    cli::cli_abort(
      c(
        "Only the following distributions are supported:",
        setNames(names(available_distributions), rep("*", length(available_distributions)))
      )
    )
  }

  if (distribution == "truncated normal") {
    cli::cli_inform(c(
      "!" = paste(
        "The Truncated Normal distribution is supported for model",
        "comparison purposes -",
        "but it is not recommended in practice due to misspecification, i.e.",
        "biased mean at low concentrations.",
        "Consider using a Gamma or Log-Normal distribution instead."
      )))
  }
  if (distribution == "normal") {
    cli::cli_inform(c(
      "!" = paste(
        "The Normal distribution is supported for model",
        "comparison purposes -",
        "but it is not recommended in practice due to misspecification, i.e.",
        "prediction of negative concentrations.",
        "Consider using a Gamma or Log-Normal distribution instead."
      )))
  }

  return(do.call(
    switch, c(distribution, as.list(c(available_distributions,-1)))
  ))
}

#' Define a prior in modeldata
#' @keywords internal
set_prior <- function(param, dist = "normal", ...) {
  prior <- c(as.list(environment()), list(...))
  prior_data <- list()
  prior_data[[paste0(param, "_prior_text")]] <- paste(
    names(prior), "=", prior,
    collapse = ", "
  )
  prior$param <- NULL
  prior$dist <- NULL
  names(prior) <- NULL
  prior_data[[paste0(param, "_prior")]] <- unlist(prior)
  return(prior_data)
}

#' Provide initialization value for a parameter based on the supplied prior with location and scale
#'
#' @description Initialization using the prior is often better than initializing
#'   with zero (and if the parameter is strictly positive, zero is not possible
#'   at all). This function provides as init value the location of the prior
#'   plus 1/4 of the scale. This ensure a positive init even if the
#'   mean is zero (useful for truncated normal priors for example.)
#'
#' @param prior Prior for parameter as provided by [set_prior()]. Should be a
#'   location and scale prior (first element location, second element scale).
#'
#' @details If the provided prior has zero variance, it is assumed that the
#'   parameter will not be sampled and an empty init is returned.
#'
#' @return Location of the prior plus 1/4 of the scale.
#' @keywords internal
init_from_location_scale_prior <- function(prior) {
  prior_data_select <- !stringr::str_detect(names(prior), pattern = "_text$")
  if (sum(prior_data_select) != 1) {
    cli::cli_warn(paste(
      "Could not init from location scale prior:",
      "Non-ambiguous prior format. Using 1e-2 as fallback."
    ), .internal = TRUE)
  }
  prior_data <- prior[prior_data_select][[1]]
  if (prior_data[2] > 0) {
    return(prior_data[1] + prior_data[2]/4)
  } else {
    return(numeric(0))
  }
}
