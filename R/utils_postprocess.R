summary_measures <- function(x) c(
  "mean" = mean(x, na.rm = TRUE),
  "median" = median(x, na.rm = TRUE),
  "sd" = sd(x, na.rm = TRUE),
  "mad" = mad(x, na.rm = TRUE),
  "lower_95" = unname(quantile(x, probs = 0.025, na.rm = TRUE)),
  "upper_95" = unname(quantile(x, probs = 0.975, na.rm = TRUE))
)

convergence_measures <- posterior::default_convergence_measures()

summarize_coefs <- function(res) {
  coef_summary <- rbindlist(list(
    alpha = as.data.table(res$fit$summary("alpha", summary_measures, convergence_measures)),
    beta = as.data.table(res$fit$summary("beta", summary_measures, convergence_measures))
  ))
  # rename beta variables with the original names
  coef_summary[, variable := c("(Intercept)", res$covariates_design)]
  return(coef_summary)
}

summarize_residuals <- function(res) {
  residuals <- as.data.table(res$fit$summary("residuals", summary_measures, convergence_measures))
  residuals <- residuals[, -c("variable")]

  # check if same length as obs_ids
  if (nrow(residuals) != nrow(res$metainfo$obs_ids)) {
    cli::cli_abort("Number of residuals does not match number of observations.")
  }

  return(cbind(res$metainfo$obs_ids, residuals))
}

#' @importFrom posterior %**%
predict_response <- function(res, X, link, interval) {

  alpha_draws <- posterior::as_draws_rvars(res$fit$draws("alpha"))$alpha
  beta_draws <- posterior::as_draws_rvars(res$fit$draws("beta"))$beta
  response <- t(alpha_draws + beta_draws %**% t(X))[,1,drop = TRUE]

  # apply link function
  if (link == "log") {
    response <- exp(response)
  }

  if (interval == "none") {
    preds <- data.table(mean = posterior::E(response))
  } else if (interval == "confidence") {
    preds <- as.data.table(posterior::summarise_draws(response, summary_measures, convergence_measures))[, -c("variable")]
  } else if (interval == "prediction") {
    cli::cli_abort("Prediction intervals are not yet supported.")
  } else {
    cli::cli_abort("Unknown interval type.")
  }
  return(preds)
}
