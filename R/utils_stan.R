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
