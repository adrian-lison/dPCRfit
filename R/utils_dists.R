#' Get shape of a Gamma distribution given its mean and sd
get_gamma_shape_alternative <- function(gamma_mean, gamma_sd) {
  gamma_shape <- (gamma_mean / gamma_sd)^2
  return(gamma_shape)
}

#' Get rate of a Gamma distribution given its mean and sd
get_gamma_rate_alternative <- function(gamma_mean, gamma_sd) {
  gamma_rate <- gamma_mean / (gamma_sd^2)
  return(gamma_rate)
}

#' Get scale of a Gamma distribution given its mean and sd
get_gamma_scale_alternative <- function(gamma_mean, gamma_sd) {
  return(1 / get_gamma_rate_alternative(gamma_mean, gamma_sd))
}

rgamma2 <- function(n, mean, cv) {
  sd = cv * mean
  shape = get_gamma_shape_alternative(mean, sd)
  rate = get_gamma_rate_alternative(mean, sd)
  return(rgamma(n = n, shape = shape, rate = rate))
}

rlnorm3 <- function(n, mean, cv) {
  sigma2 = log(1 + cv^2);
  mu = log(mean) - sigma2/2;
  return(rlnorm(n, meanlog = mu, sd = sqrt(sigma2)))
}
