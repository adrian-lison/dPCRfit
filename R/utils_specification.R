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

#' Define a normal prior in modeldata
#'
#' @param param Name of the parameter for which the prior is defined.
#' @param mu Mean.
#' @param sigma Standard deviation.
#' @param two_sigma Two times the standard deviation. Useful for defining priors
#'   via the two-sigma rule-of-thumb (approximately 95% of probability mass).
#' @param q5 Lower quantile (5%).
#' @param q95 Upper quantile (95%).
#'
#' @return Prior specification for modeldata.
#' @keywords internal
set_prior_normal <- function(param, mu = NULL, sigma = NULL, two_sigma = NULL, q5 = NULL, q95 = NULL) {
  if (!is.null(q5) && !is.null(q95)) {
    if (q5 > q95) {
      cli::cli_abort(
        "The lower quantile (q5) must be smaller than the upper quantile (q95)",
      )
    }
    mu <- (q5 + q95) / 2
    sigma <- (q95 - q5) / (2 * qnorm(0.95))
  } else {
    if (is.null(mu)) {
      cli::cli_abort(
        "mu must be supplied for the normal prior",
        .internal = TRUE
      )
    }
    if (is.null(sigma)) {
      if (is.null(two_sigma)) {
        cli::cli_abort(
          "Either sigma or two_sigma must be supplied for the normal prior",
          .internal = TRUE
        )
      } else {
        sigma <- two_sigma / 2
      }
    }
  }
  return(set_prior(param = param, dist = "normal", mu = mu, sigma = sigma))
}

#' Define a truncated normal prior in modeldata
#'
#' @description This function defines a truncated normal prior (truncation at
#'   zero) for a parameter in modeldata.
#'
#' @param param Name of the parameter for which the prior is defined.
#' @param mu Mean (not accounting for truncation).
#' @param sigma Standard deviation (not accounting for truncation).
#' @param two_sigma Two times the standard deviation (not accounting for
#'   truncation). Useful for defining priors via the two-sigma rule-of-thumb
#'   (approximately 95% of probability mass).
#' @param q5 Lower quantile (5%).
#' @param q95 Upper quantile (95%).
#'
#' @return Prior specification for modeldata.
#' @keywords internal
set_prior_trunc_normal <- function(param, mu = NULL, sigma = NULL, two_sigma = NULL, q5 = NULL, q95 = NULL) {
  if (!is.null(q5) && !is.null(q95)) {
    if (q5 > q95) {
      cli::cli_abort(
        "The lower quantile (q5) must be smaller than the upper quantile (q95)",
      )
    }
    params <- find_mu_sigma_truncnorm(q5_target = q5, q95_target = q95)
    mu <- params$mu
    sigma <- params$sigma
  } else {
    if (is.null(mu)) {
      cli::cli_abort(
        "mu must be supplied for the normal prior",
        .internal = TRUE
      )
    }
    if (is.null(sigma)) {
      if (is.null(two_sigma)) {
        cli::cli_abort(
          "Either sigma or two_sigma must be supplied for the normal prior",
          .internal = TRUE
        )
      } else {
        sigma <- two_sigma / 2
      }
    }
  }
  return(set_prior(param = param, dist = "truncated-normal", mu = mu, sigma = sigma))
}

#' Find mu and sigma of a truncated normal distribution for given quantiles
#'
#' @description This function estimates the parameters `mu` and `sigma` of a
#'   normal distribution truncated at zero, such that the 5% and 95% quantiles
#'   of the distribution match the specified target quantiles.
#'
#' @param q5_target The target quantile at 5% (lower bound).
#' @param q95_target The target quantile at 95% (upper bound).
#'
#' @return A list containing the estimated `mu` and `sigma` for the truncated
#'   normal distribution that matches the specified quantiles.
#' @keywords internal
find_mu_sigma_truncnorm <- function(q5_target, q95_target) {
  # Check if the target quantiles are valid
  if (q5_target > q95_target) {
    cli::cli_abort(paste0(
      "Provided 5% and 95% quantiles are invalid: ",
      q5_target, " - ", q95_target, "."
    ))
  }

  if (q5_target == q95_target) {
    return(list(mu = q5_target, sigma = 0))  # fixed value
  }

  # Initial guesses
  mu_start <- (q5_target + q95_target) / 2
  sigma_start <- (q95_target - q5_target) / (qnorm(0.95) - qnorm(0.05))

  # Objective function: squared error between target and actual quantiles
  objective_fn <- function(params, q5_target, q95_target) {
    mu <- params[1]
    sigma <- params[2]
    if (sigma <= 0) return(1e10)  # penalize invalid sigma

    # Compute the 5% and 95% quantiles of the truncated normal
    q5_model <- extraDistr::qtnorm(0.05, a = 0, b = Inf, mean = mu, sd = sigma)
    q95_model <- extraDistr::qtnorm(0.95, a = 0, b = Inf, mean = mu, sd = sigma)

    # Return sum of squared errors
    sum((q5_model - q5_target)^2 + (q95_model - q95_target)^2)
  }

  # Optimize to minimize the quantile error
  result <- optim(
    par = c(mu_start, sigma_start),
    fn = objective_fn,
    q5_target = q5_target,
    q95_target = q95_target,
    method = "L-BFGS-B",
    lower = c(-Inf, 1e-6),  # sigma must be positive
    control = list(fnscale = 1)
  )

  mu <- result$par[1]
  sigma <- result$par[2]

  # check if target quantiles match realized quantiles (tolerance of 1e-2)
  if (abs(extraDistr::qtnorm(0.05, mu, sigma, a = 0) - q5_target) > 1e-3 ||
      abs(extraDistr::qtnorm(0.95, mu, sigma, a = 0) - q95_target) > 1e-3) {
    cli::cli_warn(paste0(
      "90% interval of truncated normal prior could not be exactly calibrated ",
      "(specified interval: ", q5_target, " - ", q95_target, ", ",
      "realized interval: ",
      round(extraDistr::qtnorm(0.05, mu, sigma, a = 0),3),
      " - ",
      round(extraDistr::qtnorm(0.95, mu, sigma, a = 0),3),
      "). ",
      "This typically happens when the specified quantiles are rather extreme."
    ))
  }

  return(list(mu = mu, sigma = sigma))
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
init_from_location_scale_prior <- function(prior, enforce_positive = FALSE) {
  prior_data_select <- !stringr::str_detect(names(prior), pattern = "_text$")
  if (sum(prior_data_select) != 1) {
    cli::cli_warn(paste(
      "Could not init from location scale prior:",
      "Non-ambiguous prior format. Using 1e-2 as fallback."
    ), .internal = TRUE)
  }
  prior_data <- prior[prior_data_select][[1]]
  if (prior_data[2] > 0) {
    init <- prior_data[1] + prior_data[2]/4
    if (enforce_positive && init < 0) {
      init <- prior_data[2]/4
    }
    return(init)
  } else {
    return(numeric(0))
  }
}
