dPCRfit <- function(formula, data, link = "identity",
dPCRfit <- function(formula, data, link = c("identity", "log"),
                    measurements = concentration_measurements(),
                    noise = noise_dPCR(),
                    nondetect = nondetect_dPCR(),
                    prior_intercept = c(0,1),
                    prior_coefficients = c(0,1),
                    fit_opts = set_fit_opts()
                    ) {
  md <- modeldata_init()
  md$.inputs$df <- data
  link <- match.arg(link)
  md <- md +
    measurements +
    noise +
    nondetect +
    linear_regression(
      formula = formula,
      link = link,
      alpha_prior = prior_intercept,
      beta_prior = prior_coefficients
    )

  # prepare model data
  inits <- md$.init
  metainfo <- md$.metainfo
  md <- suppressWarnings(rlang::flatten(md[!(names(md) %in% c(
    ".metainfo", ".checks", ".str", ".init",
    ".sewer_data", ".sewer_assumptions"
  ))]))
  md <- md[
    stringr::str_detect(names(md), c("_prior_text"), negate = TRUE)
  ]

  # fit model
  arguments <- c(
    list(data = md),
    init = function() inits,
    fit_opts$sampler
  )

  fitting_successful <- FALSE
  result <- list(
    formula = formula,
    link = link,
    nrow = nrow(data),
    variables = colnames(md$X),
    X = md$X,
    metainfo = metainfo,
    fit_opts = fit_opts
  )

  stanmodel_instance <- cmdstanr::cmdstan_model(
    system.file(
      "stan",
      "dPCRglm.stan",
      package = "dPCRfit"
    ),
    threads = fit_opts$model$threads,
    force_recompile = fit_opts$model$force_recompile
  )

  fit_res <- tryCatch(
    {
      fit_res <- withWarnings(suppress_messages_warnings(
        do.call(stanmodel_instance$sample, arguments),
        c(
          "Registered S3 method overwritten by 'data.table'",
          "Cannot parse stat file, cannot read file: No such file or directory",
          "cannot open file '/proc/stat': No such file or directory"
        )
      ))
      if (length(fit_res$warnings) == 0) {
        fitting_successful <- TRUE
        fit_res <- fit_res$value
        # ensure that data is read in
        fit_res$draws("alpha")
      } else {
        cat("\n")
        cli::cli_warn(
          paste(
            "There was an error while fitting the model.",
            "Only the model input is returned."
          )
        )
        fit_res <- list(
          errors = unlist(lapply(
            fit_res$warnings, function(x) stringr::str_remove(x$message, "\n")
          )),
          sampler_output = fit_res$value$output()
        )
      }
      fit_res
    },
    error = function(err) {
      cat("\n")
      cli::cli_warn(c(
        paste(
          "There was an error while fitting the model.",
          "Only the model input is returned."
        ),
        err$message
      ))
      return(list(errors = err, sampler_output = NULL))
    }
  )

  result$fit <- fit_res

  if (fitting_successful) {
    tryCatch(result$coef_summary <- summarize_coefs(result),
             error = function(e) cli::cli_warn(paste(
               "Could not summarize results:", e$message
               )))
    tryCatch(result$residuals_summary <- summarize_residuals(result),
             error = function(e) cli::cli_warn(paste(
               "Could not summarize residuals:", e$message
             )))
    tryCatch(result$diagnostic_summary <- fit_res$diagnostic_summary(),
             error = function(e) cli::cli_warn(paste(
               "Could not obtain model diagnostics:", e$message
             )))
  }

  class(result) <- "dPCRfit_result"

  return(result)
}

#' @export
summary.dPCRfit_result <- function(object, ...) {
  cat("Call:\n")
  cat("dPCRfit(formula = ", format(object$formula), ", link = ", object$link, ")\n\n", sep = "")
  cat("Number of observations:", object$nrow, "\n")
  cat("\n")
  cat("Coefficients:\n")
  print(object$coef_summary)
  cat("\n")
  cat("Fitted via MCMC using ", object$fit_opts$sampler$chains, " chains with each:\n", sep = "")
  cat(object$fit_opts$sampler$iter_warmup, " warm-up iterations", "\n", sep = "")
  cat(object$fit_opts$sampler$iter_sampling, " sampling iterations", "\n", sep = "")
  cat("\n")
  cat("Diagnostics:\n")
  num_div <- sum(object$diagnostic_summary$num_divergent)
  if (num_div>0) {
    cat(paste0("Divergent transitions: ", num_div, "\n"))
  }
  num_max_treedepth <- sum(object$diagnostic_summary$num_max_treedepth)
  if (num_max_treedepth>0) {
    cat(paste0("Maximum tree depth exceeded: ", num_max_treedepth, "\n"))
  }
  num_ebfmi <- sum(object$diagnostic_summary$num_ebfmi<0.2)
  if (num_ebfmi>0) {
    cat(paste0("Low EBFMI: ", num_ebfmi, " chains\n"))
  }
  if (num_div==0 & num_max_treedepth==0 & num_ebfmi==0) {
    cat("No problems detected.\n")
  }
}

#' @export
print.dPCRfit_result <- function(object, ...) {
  cat("Call:\n")
  cat("dPCRfit(formula = ", format(object$formula), ", link = ", object$link, ")\n\n", sep = "")
  cat("Coefficients:\n")
  print(object$coef_summary[,c("variable", "mean", "median", "sd", "q5", "q95")])
}
