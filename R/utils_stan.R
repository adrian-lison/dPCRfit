#' Configure the model fitting
#'
#' @description Sets model fitting options in `dPCRfit`.
#'
#' @param sampler Which type of sampler should be used for model fitting?
#'   Currently, only [sampler_stan_mcmc()] is supported.
#' @param model Details about the model file to be used, see
#'   [model_stan_opts()].
#'
#' @return A `list` with the model fitting options.
#' @export
set_fit_opts <- function(sampler = sampler_stan_mcmc(),
                         model = model_stan_opts()) {
  opts <- as.list(environment())
  return(opts)
}

#' Use the stan MCMC sampler
#'
#' @description This option will use stan's mcmc sampler via [cmdstanr] for
#'   Markov Chain Monte Carlo (MCMC) sampling.
#'
#' @param ... Further arguments to pass to [cmdstanr].
#' @inheritParams cmdstanr::sample
#'
#' @return A `list` with settings for the MCMC sampler.
#' @export
sampler_stan_mcmc <- function(
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    max_treedepth = 15,
    step_size = 0.01,
    parallel_chains = NULL,
    threads_per_chain = 1,
    seed = 0,
    refresh = 200,
    show_messages = TRUE,
    show_exceptions = FALSE,
    ...) {
  opts <- c(as.list(environment()), list(...))
  if (opts$threads_per_chain == 1) {
    opts$threads_per_chain <- NULL
  }
  class(opts) <- "mcmc"
  return(opts)
}

#' Use stan's pathfinder variational inference algorithm
#'
#' @description This option will use stan's pathfinder algorithm via [cmdstanr]
#'   for variational inference. Sampling is very fast but potentially less
#'   exact.
#'
#' @param ... Further arguments to pass to [cmdstanr].
#' @inheritParams cmdstanr::pathfinder
#'
#' @return A `list` with settings for the pathfinder sampler.
#' @export
sampler_stan_pathfinder <- function(
    draws = 4000,
    num_paths = 4,
    single_path_draws = NULL,
    max_lbfgs_iters = NULL,
    num_elbo_draws = NULL,
    psis_resample = NULL,
    seed = 0,
    refresh = 2000,
    show_messages = TRUE,
    show_exceptions = FALSE,
    ...
  ) {
  opts <- c(as.list(environment()), list(...))
  class(opts) <- "pathfinder"
  return(opts)
}

#' Specify details of the stan model
#'
#' @param threads Should multihreading be enabled? Default is `FALSE`, as
#'   `dPCRfit` currently does not support within-chain parallelism.
#' @param force_recompile Should recompilation be forced before model fitting?
#'   If `FALSE` (default), the model is only recompiled when changes to the
#'   model code are detected. However, as the change detection is not fully
#'   reliable, it is sometimes necessary to force recompilation after having
#'   made changes to the stan code.
#'
#' @return A list with details of the stan model
#' @export
model_stan_opts <- function(threads = FALSE, force_recompile = FALSE) {
  opts <- as.list(environment())
  return(opts)
}

#' @keywords internal
get_pathfinder_inits <- function(stanmodel_instance, md, inits) {
  arguments <- c(
    list(data = md),
    init = function() inits,
    list(num_paths = 4, draws = 100, seed = 0)
  )
  return(fit_stan(
    stanmodel_instance, arguments, fit_method = "pathfinder", silent = TRUE
  ))
}

#' @keywords internal
fit_stan <- function(stanmodel_instance, arguments, fit_method, silent = FALSE) {
  if (silent) {
    arguments$show_messages <- FALSE
    arguments$show_exceptions <- FALSE
    sink(tempfile(), type = "out")
    on.exit(sink())
  }
  return(tryCatch(
    {
      fit_res <- withWarnings(suppress_messages_warnings(
        {
          if (fit_method == "mcmc") {
            do.call(stanmodel_instance$sample, arguments)
          } else if (fit_method == "pathfinder") {
            do.call(stanmodel_instance$pathfinder, arguments)
          } else {
            cli::cli_abort("Unknown sampler type")
          }
        },
        c(
          "Registered S3 method overwritten by 'data.table'",
          "Cannot parse stat file, cannot read file: No such file or directory",
          "cannot open file '/proc/stat': No such file or directory"
        )
      ))
      if (length(fit_res$warnings) == 0) {
        fit_res <- fit_res$value
        # ensure that data is read in
        fit_res$draws("alpha")
      } else {
        if (!silent) {
          cat("\n")
          cli::cli_warn(
            paste(
              "There was an error while fitting the model.",
              "Only the model input is returned."
            )
          )
        }
        fit_res <- list(
          errors = unlist(lapply(
            fit_res$warnings, function(x) stringr::str_remove(x$message, "\n")
          )),
          sampler_output = invisible(force(fit_res$value$output()))
        )
      }
      return(fit_res)
    },
    error = function(err) {
      if (!silent) {
        cat("\n")
        cli::cli_warn(c(
          paste(
            "There was an error while fitting the model.",
            "Only the model input is returned."
          ),
          err$message
        ))
      }
      return(list(errors = err, sampler_output = NULL))
    }
  ))
}

#' Compile the dPCRfit model
#'
#' @description The stan model used by dPCRfit needs to be compiled for your
#'   device. This is only necessary once after installing or updating the
#'   package.
#'
#' @param force_recompile If TRUE, the model will be recompiled even if it
#'   has already been successfully compiled.
#' @param verbose If TRUE, warnings and detailed errors from the compilation are
#'   printed. This can help to diagnose compilation issues.
#'
#' @details If the model is not successfully compiled, please
#'   ensure that `cmdstan` is properly set up and try updating it to a newer
#'   version using [cmdstanr::install_cmdstan()]. If the problem persists,
#'   please run [dPCRfit_compile(verbose = TRUE)] and post the output in
#'   a new issue on GitHub, along with your [cmdstanr::cmdstan_version()].
#'
#' @export
dPCRfit_compile <- function(force_recompile = FALSE, verbose = FALSE) {
  get_model <- function() {
    cmdstanr::cmdstan_model(
    system.file(
      "stan",
      "dPCRglm.stan",
      package = "dPCRfit"
    ),
    force_recompile = force_recompile
  )
    return(invisible(TRUE))
  }

  success <- tryCatch(
    {
      if (verbose == FALSE) {
        suppressMessages(
          get_model()
        )
      } else {
        get_model()
      }
      TRUE
    },
    error = function(e) { FALSE }
  )

  if (success) {
    cli::cli_alert("dPCRfit model compiled successfully.")
  } else {
    cli::cli_warn("The dPCRfit stan model could not be compiled.")
  }
}
