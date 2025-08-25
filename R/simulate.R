#' Simulate concentration measurements based on the dPCR partition count model
#'
#' @param n_samples How many measurements of the sample do you want to simulate?
#' @param lambda Expected concentration in the sample
#' @param m_partitions Number of partitions (vector of length n_replicates or
#'   scalar)
#' @param conversion Conversion factor, i.e. volume v * scaling factor s
#' @param prePCR_cv Pre-PCR coefficient of variation
#' @param n_replicates Number of replicates
#' @param distribution Pre-PCR variation distribution, either "gamma" or "lnorm"
#' @param get_partitions If TRUE, return a data frame with concentration,
#'   positive partitions and total partitions
#'
#' @return A vector of simulated measurements
sim_concentrations_dPCR <- function(
    n_samples = 1, lambda, m_partitions, conversion,
    prePCR_cv, n_replicates = 1, distribution = "gamma", get_partitions = FALSE) {
    # Ensure m_partitions is a vector of length n_replicates
    if (is.list(m_partitions)) {
      stopifnot(length(m_partitions[[1]]) == n_replicates)
      stopifnot(length(m_partitions) == n_samples)
    } else {
      if (length(m_partitions) == 1) m_partitions <- rep(m_partitions, n_replicates)
    }
    if (prePCR_cv==0) {
      lam <- rep(lambda, n_samples)
    } else {
      if (distribution == "gamma") {
        lam <- rgamma2(n_samples, lambda, cv = prePCR_cv)
      } else if (distribution == "lnorm") {
        lam <- rlnorm3(n_samples, lambda, cv = prePCR_cv)
      } else {
        stop("Distribution not supported.")
      }
    }
    p <- 1 - exp(-lam * conversion)
    sim_res <- mapply(function(p, m_partitions) {
      sim_drop <- mapply(function(mp) rbinom(1, mp, p), m_partitions)
      sim_positive_partitions <- sum(sim_drop)
      sim_total_partitions <- sum(m_partitions)
      sim_conc <- - log(1 - (sim_positive_partitions / sim_total_partitions)) * (1/conversion)
      return(list(
        concentration = sim_conc,
        positive_partitions = sim_positive_partitions/n_replicates, # average over replicates
        total_partitions = sim_total_partitions/n_replicates # average over replicates
      ))
    }, p = p, m_partitions = m_partitions, SIMPLIFY = FALSE)
    if (get_partitions) {
      return(data.frame(
        concentration = sapply(sim_res, function(x) x$concentration),
        positive_partitions = sapply(sim_res, function(x) x$positive_partitions),
        total_partitions = sapply(sim_res, function(x) x$total_partitions)
      ))
    } else {
      return(sapply(sim_res, function(x) x$concentration))
    }
  }

#' Simulate number of valid partitions across dPCR runs
#'
#' @param n_samples Number of different samples
#' @param n_replicates Number of replicates per sample
#' @param max_partitions Maximum number of partitions
#' @param partition_loss_logit_mean Mean proportion of invalid partitions
#'   (logit-level)
#' @param partition_loss_logit_sd Sd of proportion of invalid partitions
#'   (logit-level)
#' @param partition_loss_max Upper bound for partition loss proportion
#'
#' @return Simulated partition numbers, structured as a list with one element
#'   per sample, where each element is a vector giving the partition numbers
#'   across replicates
sim_partitions <- function(n_samples, n_replicates, max_partitions, partition_loss_logit_mean, partition_loss_logit_sd, partition_loss_max) {
  lapply(1:n_samples, function(x) extraDistr::rtbinom(n_replicates, size = max_partitions, prob = 1-plogis(rnorm(n_replicates, mean = partition_loss_logit_mean, sd = partition_loss_logit_sd)), a = max_partitions * (1-partition_loss_max)))
}
