#' Fit a dPCR-specific regression model
#'
#' @description Fits regression models to gene concentration measurements from
#'   digital PCR (dPCR) using a dPCR-specific likelihood that accounts for
#'   concentration-dependent error and non-detects. Supports multiple linear
#'   regression with identity or log link, incorporating prior information on
#'   lab parameters. Models are fitted via Stan using the 'cmdstanr' interface.
#'
#' @param formula A symbolic description of the model to be fitted, similar to a
#'   formula provided to `lm`. This is where you specify the covariates.
#' @param data A data.frame containing the variables in the model.
#' @param link The link function to use. Can be "identity" or "log".
#' @param measurements Assumptions about the measurement process and laboratory
#'   parameters, provided through the [concentration_measurements()] component.
#' @param noise Assumptions about the noise in the data, provided through the
#'   [noise_dPCR()] or the [noise_constant_cv()] / [noise_constant_var()]
#'   component.
#' @param nondetect Assumptions about the non-detects in the data, provided
#'   through the [nondetect_dPCR()] or the [nondetect_none()] component.
#' @param prior_intercept A vector of two values (mean, sd) specifying the
#'   normal prior for the intercept.
#' @param prior_coefficients A vector of two values (mean, sd) specifying the
#'   normal prior for the regression coefficients.
#' @param fit_opts A list of options for the fitting process. See
#'   `set_fit_opts()` for details.
#'
#' @return An object of class "dPCRfit_result" containing the fitted model and
#'   additional details.
#'
#' @import stats
#' @import data.table
#'
#' @export
dPCRfit <- function(formula, data, link = c("identity", "log"),
                    measurements = concentration_measurements(),
                    noise = noise_dPCR(),
                    nondetect = nondetect_dPCR(),
                    prior_intercept = c(0,1),
                    prior_coefficients = c(0,1),
                    fit_opts = set_fit_opts()
                    ) {
  md <- modeldata_init()
  data <- setDT(data)
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

  # run checks
  for (check in md$.checks) {
    check(md, data)
  }

  # prepare model data
  inits <- md$.init
  metainfo <- md$.metainfo
  md <- suppressWarnings(rlang::flatten(md[!(names(md) %in% c(
    ".metainfo", ".checks", ".str", ".init",
    ".sewer_data", ".sewer_assumptions", ".inputs"
  ))]))
  md <- md[
    stringr::str_detect(names(md), c("_prior_text"), negate = TRUE)
  ]

  # fit model
  result <- list(
    formula = formula,
    link = link,
    nrow = nrow(data),
    covariates = all.vars(rlang::f_rhs(formula)),
    covariates_design = colnames(md$X),
    covariates_df = data[, .SD, .SDcols = all.vars(rlang::f_rhs(formula))],
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

  arguments <- c(
    list(data = md),
    init = function() inits,
    fit_opts$sampler
  )

  # pathfinder initialization for mcmc
  if (class(fit_opts$sampler) == "mcmc") {
    tryCatch(
      {
        cat("Initializing chains via pathfinder...\n")
        pathfind_init <- get_pathfinder_inits(stanmodel_instance, md, inits)
        stopifnot(!"errors" %in% names(pathfind_init))
        options(cmdstanr_warn_inits = FALSE)
        arguments$init <- pathfind_init$draws()
      },
      error = function(e) {
        cli::cli_warn(paste(
          "Pathfinder initialization failed.",
          "Falling back to default initialization."
        ))
      }
    )
  }

  result$fit <- fit_stan(stanmodel_instance, arguments, fit_method = class(fit_opts$sampler))
  options(cmdstanr_warn_inits = TRUE)

  if (!"errors" %in% names(result$fit)) {
    tryCatch(result$coef_summary <- summarize_coefs(result),
             error = function(e) cli::cli_warn(paste(
               "Could not summarize results:", e$message
               )))
    tryCatch(result$residuals_summary <- summarize_residuals(result),
             error = function(e) cli::cli_warn(paste(
               "Could not summarize residuals:", e$message
             )))
    if (class(fit_opts$sampler) == "mcmc") {
      tryCatch(result$diagnostic_summary <- suppressMessages(result$fit$diagnostic_summary()),
               error = function(e) cli::cli_warn(paste(
                 "Could not obtain model diagnostics:", e$message
               )))
    }
  }

  class(result) <- "dPCRfit_result"

  return(result)
}

#' Provide a summary of a dPCR Model Fit
#'
#' @description Summarizes a dPCR Model Fit, including coefficient estimates and
#'   model diagnostics.
#'
#' @export
summary.dPCRfit_result <- function(object, ...) {
  cat("Call:\n")
  cat("dPCRfit(formula = ", format(object$formula), ", link = ", object$link, ")\n\n", sep = "")
  cat("Number of observations:", object$nrow, "\n")
  cat("\n")

  if ("errors" %in% names(object$fit)) {
    cat(paste(
      "The model was not successfully fitted,",
      "see $fit$errors and $fit$sampler_output for details."
    ))
  } else {
    cat("Coefficients:\n")
    # print object$coef_summary as character table, not notebook output
    cat(format_table(object$coef_summary))
    cat("\n")
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
}

#' Print method for dPCR Model Fits
#'
#' @description Prints a short summary of a dPCR Model Fit
#'
#' @export
print.dPCRfit_result <- function(object, ...) {
  cat("Call:\n")
  cat("dPCRfit(formula = ", format(object$formula), ", link = ", object$link, ")\n\n", sep = "")

  if ("errors" %in% names(object$fit)) {
    cat(paste(
      "The model was not successfully fitted,",
      "see $fit$errors and $fit$sampler_output for details."
      ))
  } else {
    cat("Coefficients:\n")
    cat(format_table(object$coef_summary[,c("variable", "mean", "median", "sd", "q5", "q95")]))
  }
}

#' Predict method for dPCR Model Fits
#'
#' @description Predict values based on a fitted dPCR-specific regression model.
#'
#' @param object Object of class "dPCRfit_result"
#' @param newdata An optional data frame in which to look for covariates with
#'   which to predict. If omitted, the fitted values are used.
#' @param interval Type of interval calculation. Can be abbreviated. Currently,
#'   only "none" and "confidence" are supported. "None" will only give the mean
#'   response. Confidence will give various summary statistics of the response.
#'   Prediction intervals (i.e. including dPCR measurement noise) will be
#'   added in the future.
#' @param keep_data Logical. If TRUE, the newdata is returned together with the
#'   predictions.
#'
#' @return A data.table with the predicted values. If keep_data is TRUE, the
#'   covariates are included in this data.table.
#'
#' @export
predict.dPCRfit_result <- function(object, newdata, interval = c("none", "confidence", "prediction"), keep_data = FALSE) {

  if ("errors" %in% names(object$fit)) {
    cli::cli_abort(paste(
      "The model was not successfully fitted,",
      "see $fit$errors and $fit$sampler_output for details."
    ))
  }

  if (missing(newdata) || is.null(newdata)) {
    newdata <- data.table::copy(object$covariates_df)
  } else if (!(is.data.frame(newdata) || is.matrix(newdata))) {
    cli::cli_abort("newdata must be a data.frame or matrix")
  } else if (!all(object$covariates %in% colnames(newdata))) {
    cli::cli_abort(paste(
      "newdata must contain the following covariates:",
      paste(object$covariates, collapse = ", ")
      ))
  }
  setDT(newdata)
  newdata[[all.vars(rlang::f_lhs(object$formula))]] <- 0 # dummy response
  newdata <- newdata[, .SD, .SDcols = all.vars(object$formula)]
  if (any(duplicated(newdata))) {
    cli::cli_warn("Duplicated rows in newdata are removed.")
    newdata <- newdata[!duplicated(newdata), ]
  }
  for (col in object$covariates) {
    col_data <- object$covariates_df[[col]]
    if (is.factor(col_data) || is.character(col_data)) {
      if (is.character(col_data)) {
        col_data <- factor(col_data)
      }
      if (!all(unique(newdata[[col]]) %in% levels(col_data))) {
        cli::cli_abort(paste(
          "Variable", col,
          "in newdata contains level(s) not seen by the model:",
          paste(setdiff(unique(newdata[[col]]), levels(col_data)), collapse = ", ")
          ))
      }
      newdata[[col]] <- factor(
        newdata[[col]],
        levels = levels(col_data),
        ordered = is.ordered(col_data)
        )
    }
  }

  X <- model.matrix(object$formula, data = newdata)
  # if first column of X corresponds to an intercept, remove it
  if (all(X[, 1] == 1)) {
    X <- subset(X, select = -1)
  }
  newdata[[all.vars(rlang::f_lhs(object$formula))]] <- NULL # remove dummy response
  colnames(X) <- object$covariates_design

  interval <- match.arg(interval)
  preds <- predict_response(object, X, object$link, interval)

  if (keep_data) {
    preds <- cbind(newdata, preds)
  }

  return(preds)
}
