#' Model concentration measurements
#'
#' @param measurements A `data.frame` with each row representing one dPCR
#'   measurement. Must have at least a column with sample ids and a column with
#'   concentration measurements.
#' @param composite_window Over how many days has each measured sample been
#'   collected? If 1 (default), samples represent single days. If larger than 1,
#'   samples are assumed to be equivolumetric composites over several dates. In
#'   this case, the supplied dates represent the last day included in each
#'   sample.
#' @param id_col Name of the column containing the id of the sample.
#' @param concentration_col Name of the column containing the measured
#'   concentrations.
#' @param replicate_col Name of the column containing the replicate ID of each
#'   measurement. This is used to identify multiple measurements made of a
#'   sample from the same date. Should be `NULL` if only one measurement per
#'   sample was made.
#' @param n_averaged The number of replicates over which the measurements have
#'   been averaged. This is typically used as an alternative to providing
#'   several replicates per sample (i.e. the concentration provided in the
#'   `measurements` `data.frame` is the average of several replicates).
#'   Can be either a single number (it is then assumed that the number of
#'   averaged replicates is the same for each observation) or a vector (one
#'   value for each observation).
#' @param n_averaged_col Name of the column in the `measurements` data.frame
#'   containing the number of replicates over which the measurements have been
#'   averaged. This is an alternative to specifying `n_averaged`.
#' @param total_partitions_col Name of the column in the `measurements`
#'   data.frame containing the total number of partitions (e.g. droplets for
#'   ddPCR) in the dPCR reaction of each measurement. Only applies to
#'   concentration measurements obtain via dPCR. Can be used by the
#'   [noise_estimate_dPCR()] and [LOD_estimate_dPCR()] modeling components. Note
#'   that this is really the
#'   *total* number of partitions, not just the number of positive partitions.
#' @param distribution Parametric distribution for concentration measurements.
#'   Currently supported are "gamma" (default and recommended), "log-normal",
#'   "truncated normal", and "normal". The "truncated normal" and "normal"
#'   options are not recommended for use in practice.
#'
#' @export
concentration_measurements <-
  function(measurements = NULL,
           id_col = "sample_id",
           concentration_col = "concentration",
           replicate_col = NULL,
           n_averaged = 1,
           n_averaged_col = NULL,
           total_partitions_col = NULL,
           distribution = "gamma") {
    model_component("concentration_measurements", {
      if (is.null(measurements)) {
        if (is.null(modeldata$.inputs$df)) {
          cli::cli_abort("Please provide a `data.frame` with measurements.")
        } else {
          measurements <- modeldata$.inputs$df
        }
      }

      measurements <- select_required_cols(
        measurements,
        required_cols = list(
          id_col, concentration_col, replicate_col,
          n_averaged_col, total_partitions_col
        ),
        col_names = c(
          "sample_id", "concentration", "replicate_id",
          "n_averaged", "total_partitions"
        )
      )

      measurements <- measurements[!is.na(concentration) & !is.na(sample_id), ]
      measurements[, concentration := as.numeric(concentration)]

      modeldata$.inputs$id_col <- id_col

      if (nrow(measurements)==0) {
        cli::cli_abort(
          "The provided measurements `data.frame` contains no valid measurements."
        )
      }

      # measurements and samples
      modeldata$n_measured <- nrow(measurements)
      modeldata$n_samples <- length(unique(measurements[["sample_id"]]))
      modeldata$.metainfo$sample_ids <- sort(unique(measurements[["sample_id"]]))
      modeldata$measure_to_sample <- sapply(measurements[["sample_id"]], function(x) {
        which(x == modeldata$.metainfo$sample_ids)[[1]]
      })
      modeldata$measured_concentrations <- measurements[["concentration"]]

      # explicit replicates
      if (!is.null(replicate_col)) {
        modeldata$.metainfo$obs_ids <- measurements[, c("sample_id", "replicate_id")]
        modeldata$replicate_ids <- as.integer(measurements[["replicate_id"]])
      } else {
        modeldata$.metainfo$obs_ids <- measurements[, c("sample_id")]
      }

      # number of averaged technical replicates
      if (!is.null(n_averaged_col)){
        modeldata$n_averaged <- as.numeric(measurements[["n_averaged"]])
        if (any(is.na(modeldata$n_averaged))) {
          cli::cli_abort(paste0(
            "The column `", n_averaged_col, "` contains missing ",
            "values for some of the observed measurements."
          ))
        }
      } else if (length(n_averaged) == 1) {
        modeldata$n_averaged <- rep(n_averaged, modeldata$n_measured)
      } else if (length(n_averaged) == modeldata$n_measured) {
        modeldata$n_averaged <- n_averaged
      } else {
        cli::cli_abort(paste(
          "The length of `n_averaged` must be either 1 or equal to the",
          "number of samples."
        ))
      }

      # total valid partitions in PCR run
      if (!is.null(total_partitions_col)) {
        modeldata$dPCR_total_partitions <- as.integer(
          measurements[["total_partitions"]]
        )
      } else {
        modeldata$dPCR_total_partitions <- rep(0, modeldata$n_measured)
      }

      # distribution of non-zero measurements
      modeldata$obs_dist = set_distribution(distribution)

      return(modeldata)
    })
  }

#' Model measurement noise (internal helper function)
#'
#' @description This helper function is called from specific noise modeling
#'   function. [noise_constant_cv()] is a constant coefficient of variation
#'   model, [noise_constant_var()] is a constant variance model, and
#'   [noise_dPCR()] is a noise model specialized for digital PCR (`cv_type =
#'   "dPCR"`).
#'
#' @param cv_prior_mu Prior (mean) on the coefficient of variation of
#'   concentration measurements.
#' @param cv_prior_sigma Prior (standard deviation) on the coefficient of
#'   variation of concentration measurements.
#' @param cv_type One out of "constant" (default), "constant_var", or "dPCR". If
#'   "constant", the coefficient of variation is estimated as a constant/single
#'   parameter for all observations. If "dPCR", the coefficient of variation is
#'   modeled as a function of the expected concentration according to the
#'   statistical properties of dPCR. In particular, this model predicts a higher
#'   coefficient of variation at smaller concentrations, which often leads to a
#'   better model fit. If "constant_var", not the coefficient of variation but
#'   the variance of measurements is modeled as constant. This is usually a
#'   misspecification and is only supported for comparison purposes.
#' @param total_partitions_prior_mu Prior (mean) on the total number of
#'   partitions in the dPCR reaction.
#' @param total_partitions_prior_sigma Prior (standard deviation) on the total
#'   number of partitions in the dPCR reaction. If this is set to zero, the
#'   total number of partitions will be fixed to the prior mean and not
#'   estimated.
#' @param total_partitions_observe If TRUE, the total number of partitions is
#'   taken from the supplied measurements `data.frame`. This requires that the
#'   argument `total_partitions_col` is specified in [concentrations_observe()].
#' @param partition_variation_prior_mu Prior (mean) on the coefficient of
#'   variation of the total number of partitions in the dPCR reaction. Usually,
#'   the maximum number of partitions possible for a given dPCR chip is not
#'   reached, i.e. a certain number of partitions is lost. This loss varies
#'   between PCR runs, and is modeled as log-normal distributed.
#' @param partition_variation_prior_sigma Prior (standard deviation) on the
#'   coefficient of variation of the total number of partitions in the dPCR
#'   reaction. If this is set to zero, the partition variation will be fixed to
#'   the prior mean and not estimated.
#' @param volume_scaled_prior_mu Prior (mean) on the conversion factor
#'   (partition volume multiplied with the scaling of concentration in the
#'   assay) for the dPCR reaction. See details for further explanation.
#' @param volume_scaled_prior_sigma Prior (standard deviation) on the conversion
#'   factor (partition volume multiplied with the scaling of concentration in
#'   the assay) for the dPCR reaction. If this is set to zero, the conversion
#'   factor will be fixed to the prior mean and not estimated.
#' @param prePCR_noise_type The parametric distribution to assume for noise
#'   before the PCR assay. Currently supported are "log-normal" and "gamma". The
#'   choice of the parametric distribution typically makes no relevant
#'   difference for the noise model, but can make a relevant difference for the
#'   LOD model if [LOD_estimate_dPCR()] is used.
#'
#' @param use_taylor_approx If TRUE (default), a Taylor expansion approximation
#'   is used to estimate the CV of measurements under pre-PCR noise. The
#'   approximation is very accurate, unless concentrations are extremely high
#'   (so high that the quality of the measurements from dPCR would anyway be
#'   questionable).
#'
#' @keywords internal
noise_ <-
  function(cv_prior_mu = 0,
           cv_prior_sigma = 1,
           cv_type = "constant_cv",
           total_partitions_prior_mu = NULL,
           total_partitions_prior_sigma = NULL,
           total_partitions_observe = NULL,
           partition_variation_prior_mu = NULL,
           partition_variation_prior_sigma = NULL,
           volume_scaled_prior_mu = NULL,
           volume_scaled_prior_sigma = NULL,
           prePCR_noise_type = "log-normal",
           use_taylor_approx = TRUE,
           modeldata) {
    modeldata$nu_upsilon_a_prior <- set_prior(
      "nu_upsilon_a", "truncated normal",
      mu = cv_prior_mu, sigma = cv_prior_sigma
    )
    modeldata$.init$nu_upsilon_a <- 0.1 # 10% coefficient of variation

    if (cv_type == "constant_cv") {
      modeldata$total_partitions_observe <- FALSE
      modeldata$cv_type <- 0
      modeldata$nu_upsilon_b_mu_prior <- numeric(0)
      modeldata$nu_upsilon_b_cv_prior <- numeric(0)
      modeldata$.init$nu_upsilon_b_mu <- numeric(0)
      modeldata$.init$nu_upsilon_b_cv <- numeric(0)
      modeldata$.init$nu_upsilon_b_noise_raw <- numeric(0)
      modeldata$nu_upsilon_c_prior <- numeric(0)
      modeldata$.init$nu_upsilon_c <- numeric(0)
      modeldata$cv_pre_type <- numeric(0)
      modeldata$cv_pre_approx_taylor <- numeric(0)
    } else if (cv_type == "dPCR") {
      modeldata$cv_type <- 1

      if (total_partitions_observe) {
        modeldata$total_partitions_observe <- TRUE
        modeldata$nu_upsilon_b_mu_prior <- set_prior(
          "nu_upsilon_b_mu", "dummy prior", mu = 0, sigma = 0
        )
        modeldata$nu_upsilon_b_cv_prior <- numeric(0)
        modeldata$.init$nu_upsilon_b_mu <- numeric(0)
        modeldata$.init$nu_upsilon_b_cv <- numeric(0)
        modeldata$.init$nu_upsilon_b_noise_raw <- numeric(0)
        modeldata$.checks$check_total_partitions_col <- function(md, ...) {
          if (!"dPCR_total_partitions" %in% names(md)) {
            cli::cli_abort(paste0(
              "You specified `total_partitions_observe = TRUE`, which requires ",
              "a column with the number of total partitions in the PCR for ",
              "each sample in your data. Please specify such a column via the",
              "`total_partitions_col` argument in ",
              cli_help("concentrations_observe"), "."
            ))
          }
        }
      } else {
        modeldata$total_partitions_observe <- FALSE

        modeldata$nu_upsilon_b_mu_prior <- set_prior(
          "nu_upsilon_b_mu", "truncated normal",
          mu = total_partitions_prior_mu * 1e-4, # scaling by 1e-4 for numerical reasons
          sigma = total_partitions_prior_sigma * 1e-4
        )
        modeldata$.init$nu_upsilon_b_mu <- init_from_location_scale_prior(
          modeldata$nu_upsilon_b_mu_prior
        )

        modeldata$nu_upsilon_b_cv_prior <- set_prior(
          "nu_upsilon_b_cv", "truncated normal",
          mu = partition_variation_prior_mu,
          sigma = partition_variation_prior_sigma
        )
        modeldata$.init$nu_upsilon_b_cv <- init_from_location_scale_prior(
          modeldata$nu_upsilon_b_cv_prior
        )
        modeldata$.init$nu_upsilon_b_noise_raw <- rep(0, modeldata$n_measured)
      }

      modeldata$nu_upsilon_c_prior <- set_prior(
        "nu_upsilon_c", "truncated normal",
        mu = volume_scaled_prior_mu * 1e+5, # scaling by 1e+5 for numerical reasons
        sigma = volume_scaled_prior_sigma * 1e+5
      )
      modeldata$.init$nu_upsilon_c <- init_from_location_scale_prior(
        modeldata$nu_upsilon_c_prior
      )

      if (prePCR_noise_type == "gamma") {
        modeldata$cv_pre_type <- 0
      } else if (prePCR_noise_type %in% c("log-normal", "lognormal")) {
        modeldata$cv_pre_type <- 1
      } else {
        cli::cli_abort(paste0(
          "`prePCR_noise_type = ", prePCR_noise_type, "` not supported.",
          "Available options: 'gamma', `log-normal`."
        ))
      }
      modeldata$cv_pre_approx_taylor <- use_taylor_approx

    } else if (cv_type == "constant_var") {
      modeldata$total_partitions_observe <- FALSE
      modeldata$cv_type <- 2
      modeldata$nu_upsilon_b_mu_prior <- numeric(0)
      modeldata$nu_upsilon_b_cv_prior <- numeric(0)
      modeldata$.init$nu_upsilon_b_mu <- numeric(0)
      modeldata$.init$nu_upsilon_b_cv <- numeric(0)
      modeldata$.init$nu_upsilon_b_noise_raw <- numeric(0)
      modeldata$nu_upsilon_c_prior <- numeric(0)
      modeldata$.init$nu_upsilon_c <- numeric(0)
      modeldata$cv_pre_type <- numeric(0)
      modeldata$cv_pre_approx_taylor <- numeric(0)
    } else {
      cli::cli_abort(
        paste0(
          "Noise type `", cv_type, "` not supported. Available options: ",
          "'constant', `dPCR`, `constant_var`."
        ),
      )
    }
    return(modeldata)
  }

#' Model measurement noise with constant coefficient of variation
#'
#' @description This option estimates measurement noise via a constant
#'   coefficient of variation model.
#'
#' @description For a non-constant coefficient of variation model, see
#'   [noise_dPCR()].
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements (`cv`): `Truncated normal`
#'
#' @inheritParams noise_
#' @export
noise_constant_cv <- function(cv_prior_mu = 0, cv_prior_sigma = 1) {
    model_component("noise_constant_cv", {
      noise_(
        cv_prior_mu = cv_prior_mu,
        cv_prior_sigma = cv_prior_sigma,
        cv_type = "constant_cv",
        modeldata = modeldata
      )
    })
  }

#' Model measurement noise with constant variance
#'
#' @description This option estimates measurement noise via a constant variance
#'   model.This is a misspecification for dPCR data and is only supported for
#'   comparison purposes.
#'
#' @description For a constant coefficient of variation model, see
#'   [noise_constant_cv()], and for a non-constant coefficient of
#'   variation model, see [noise_dPCR()].
#'
#' @details Note that although this model keeps the variance constant, the prior
#'   for the measurement noise is still in terms of the (average) coefficient of
#'   variation (CV). This makes prior specification easier since the CV is
#'   unitless.
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements: `Truncated normal`
#'
#' @inheritParams noise_
#' @export
noise_constant_var <-
  function(cv_prior_mu = 0,
           cv_prior_sigma = 1,
           warn = TRUE) {
    if (warn) {
      cli::cli_inform(c("!" = paste0(
        "You have specified ",
        cli_help("noise_constant_var"),
        " as the model component for measurement noise.",
        " Note that modeling a constant variance",
        " is likely a model misspecification and should only be used for",
        " comparison purposes with better models like ",
        cli_help("noise_constant_cv"), " or ",
        cli_help("noise_dPCR"), ".",
        " You can specify ",
        "{.code noise_constant_var(warn=FALSE)} to disable this warning."
      )))
    }
    model_component("noise_constant_var", {
      noise_(
        cv_prior_mu = cv_prior_mu,
        cv_prior_sigma = cv_prior_sigma,
        cv_type = "constant_var",
        modeldata = modeldata
      )
    })
  }

#' Model measurement noise for digital PCR data
#'
#' @description This option models concentration measurements using a
#'   coefficient of variation model specialized for digital PCR (e.g. digital
#'   droplet PCR). The coefficient of variation is modeled as a function of the
#'   expected concentration according to the statistical properties of dPCR.
#'
#' @param cv_prior_mu Prior (mean) on the coefficient of variation of
#'   concentration measurements. Note that in contrast to using
#'   [noise_estimate()], this does *not* include the technical noise of the
#'   digital PCR. This is because the dPCR noise is explicitly modeled (using
#'   the number of partitions and conversion factor).
#'
#' @details The conversion factor (see `volume_scaled_prior_mu`,
#'   `volume_scaled_prior_sigma`) is the partition volume v multiplied with a
#'   scaling factor s. The scaling factor accounts concentration differences
#'   between the sample and the reaction mix, for example due to extraction or
#'   adding of reagents. For example, if the partition volume is 4.5e-7 mL and the
#'   scaling factor is 100:3 (i.e. 100 gc/mL in the original
#'   sample correspond to 3 gc/mL in the PCR reaction), then the
#'   overall conversion factor is 4.5e-7 * 100 / 3 = 1.5e-5.
#'
#' @details The priors of this component have the following functional form:
#' - coefficient of variation of concentration measurements (`cv`): `Truncated normal`
#' - mean number of partitions in dPCR: `Truncated normal`
#' - coefficient of variation of number of partitions in dPCR: `Truncated normal`
#' - conversion factor for dPCR: `Truncated normal`
#'
#' @inheritParams noise_
#' @export
noise_dPCR <-
  function(cv_prior_mu = 0,
           cv_prior_sigma = 1,
           total_partitions_prior_mu = 20000,
           total_partitions_prior_sigma = 5000,
           total_partitions_observe = FALSE,
           partition_variation_prior_mu = 0,
           partition_variation_prior_sigma = 0.05,
           volume_scaled_prior_mu = 1e-5,
           volume_scaled_prior_sigma = 4e-5,
           prePCR_noise_type = "log-normal",
           use_taylor_approx = TRUE) {
    model_component("noise_dPCR", {
      noise_(
        cv_prior_mu = cv_prior_mu,
        cv_prior_sigma = cv_prior_sigma,
        cv_type = "dPCR",
        total_partitions_prior_mu = total_partitions_prior_mu,
        total_partitions_prior_sigma = total_partitions_prior_sigma,
        total_partitions_observe = total_partitions_observe,
        partition_variation_prior_mu = partition_variation_prior_mu,
        partition_variation_prior_sigma = partition_variation_prior_sigma,
        volume_scaled_prior_mu = volume_scaled_prior_mu,
        volume_scaled_prior_sigma = volume_scaled_prior_sigma,
        prePCR_noise_type = prePCR_noise_type,
        use_taylor_approx = use_taylor_approx,
        modeldata = modeldata
      )
    })
  }

#' Do not model non-detects
#'
#' @description This option drops all zero measurements from the likelihood.
#'
#' @details Dropping zero measurements implicitly conditions on positive
#'   measurements. This can lead to an upward bias when estimating low
#'   concentrations.
#'
#' @export
nondetect_none <- function() {
  model_component("nondetect_none", {
    modeldata$LOD_model <- 0
    modeldata$LOD_scale <- numeric(0)
    modeldata$LOD_drop_prob <- 0
    return(modeldata)
  })
}

#' Model non-detects in digital PCR data
#'
#' @description Concentrations below a certain threshold may not be detectable
#'   and thus erroneously measured as 0. This option adjusts for a
#'   non-detects based on the statistical properties of digital PCR (dPCR)
#'   and includes zero measurements in the likelihood.
#'
#' @description In effect, zero measurements provide a signal that the
#'   concentration in the respective sample was likely below the limit of
#'   detection, but we don't know what the exact concentration was.
#'
#' @param drop_prob Probability for non-detection below which likelihood
#'   contributions of observed concentrations are dropped. This avoids numerical
#'   issues of the non-detection model at high concentrations (very small
#'   non-detection probabilities) that can otherwise affect sampling speed.
#'   Since these likelihood contributions will be virtually zero for almost all
#'   samples anyway, parameter estimates are practically not affected.
#'
#' @details Non-detects are modeled using a hurdle model. The model
#'   uses the number of partitions in the dPCR reaction and the conversion
#'   factor as defined and estimated by [noise_estimate_dPCR()]. It can
#'   therefore only be used together with `noise = noise_estimate_dPCR()` in
#'   [concentration_measurements()].
#'
#' @export
#'
#' @family {LOD models}
nondetect_dPCR <- function(drop_prob = 1e-10) {
 model_component("nondetect_dPCR", {
    modeldata$LOD_model <- 2
    modeldata$LOD_scale <- numeric(0)

    if (is.null(modeldata$cv_type) || modeldata$cv_type != 1) {
      cli::cli_abort(paste0(
        "To use ",
        cli_help("nondetect_dPCR"), ", you must also use ",
        cli_help("noise_dPCR"), "."
      ))
    }

    modeldata$LOD_drop_prob <- drop_prob

    return(modeldata)
  })
}

