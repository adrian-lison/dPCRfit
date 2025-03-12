# Simulation of measurements based on the droplet count model
sim_concentrations_dPCR <- function(
    n_samples = 1, lambda, total_partitions, conversion,
    prePCR_cv, n_replicates = 1, distribution = "gamma") {
  # simulate pre-PCR variation
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
  concentrations <- sapply(p, function(p) {
    sim_part <- rbinom(n_replicates, total_partitions, p)
    sim_conc <- - log(1 - (sim_part / total_partitions)) * (1/conversion)
    return(mean(sim_conc))
  })
  return(concentrations)
}
