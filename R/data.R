#' Simulated dPCR data from a linear model
#'
#' This is exemplary simulated dPCR data with an underlying linear relationship
#' between the expected concentration in the sample and a covariate
#' `temperature`. The data is stored in a long format, with each row
#' representing a single measurement.
#'
#' The dPCR measurements were simulated assuming a dPCR assay with 25000
#' partitions, a conversion factor of 1.73e-5, and a pre-PCR coefficient of
#' variation of 10%.
#'
#' @format A data frame with 102 rows and 5 variables:
#' \describe{
#'   \item{sample_id}{A character vector with the sample ID.}
#'  \item{replicate_id}{A character vector with the replicate ID.}
#'   \item{biomass}{A numeric vector with the total biomass in the sampled location.}
#'   \item{lambda}{A numeric vector with the true concentration values.}
#'   \item{concentration}{A numeric vector with the measured concentrations.}
#' }
#' @source Simulated dataset
"dPCR_linear_simulated"
