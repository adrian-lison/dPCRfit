% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_specification.R
\name{init_from_location_scale_prior}
\alias{init_from_location_scale_prior}
\title{Provide initialization value for a parameter based on the supplied prior with location and scale}
\usage{
init_from_location_scale_prior(prior)
}
\arguments{
\item{prior}{Prior for parameter as provided by \code{\link[=set_prior]{set_prior()}}. Should be a
location and scale prior (first element location, second element scale).}
}
\value{
Location of the prior plus 1/4 of the scale.
}
\description{
Initialization using the prior is often better than initializing
with zero (and if the parameter is strictly positive, zero is not possible
at all). This function provides as init value the location of the prior
plus 1/4 of the scale. This ensure a positive init even if the
mean is zero (useful for truncated normal priors for example.)
}
\details{
If the provided prior has zero variance, it is assumed that the
parameter will not be sampled and an empty init is returned.
}
\keyword{internal}
