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
#' @description This option will use stan's NUTS sampler via [cmdstanr] for
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
