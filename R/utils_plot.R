#' Plot prior distributions for model components
#'
#' @description Plotting prior distributions as specified in the different model
#'   components can help to identify accidental prior misspecification and to
#'   understand the implications of specific prior choices.
#'
#' @param component A model component object. Currently the following components
#' are supported:
#'  - [noise_dPCR()]: dPCR noise model component
#' @param type Type of prior to plot. Currently the following priors are
#'   supported:
#'  - component [noise_dPCR()]: "partitions"
#' @param n_draws Number of draws to simulate from the prior distributions.
#' @param show_draws Number of example draws to show in the plots.
#' @param seed Seed for random number generation to ensure reproducibility.
#'
#' @return A ggplot object showing the prior distributions for the specified
#'   component and type.
#' @export
#' @import ggplot2
#' @import cowplot
#' @examples
#' noise_model <- noise_dPCR()
#' plot_prior(noise_model, type = "partitions")
plot_prior <- function(component, type, n_draws = 1000, show_draws = 50, seed = 0) {
 if (class(component) != "model_component") {
   cli::cli_abort("Please provide a model_component object.")
  }

  supported_components <- list(
    noise_dPCR = c("partitions")
  )

  if (!component$name %in% names(supported_components)) {
    cli::cli_abort(c(
      "Component {.val {component$name}} currently not supported for plotting of priors.",
      "Supported components are: {.val {names(supported_components)}}"
    ))
  }

  if(!type %in% supported_components[[component$name]]) {
    cli::cli_abort(c(
      "Prior of type {.val {type}} currently not supported for component {.val {component$name}}.",
      "Supported types are: {.val {supported_components[[component$name]]}}"
    ))
  }

  if (component$name == "noise_dPCR" && type == "partitions") {
    if (component$args$partitions_observe == TRUE) {
      cli::cli_abort(c(
        paste("Component {.val {component$name}} does not support plotting of",
              "the prior for {.val {type}} when partitions are observed."),
        paste("To plot a prior for the number of partitions, please set",
              "{.code partitions_observe = FALSE} in the component arguments.")
      ))
    }
    pior_md <- modeldata_init() + component
    prior_plot <- plot_prior_partitions(
      max_partitions_prior = pior_md$max_partitions_prior$max_partitions_prior * 1e4,
      partition_loss_mu_prior = pior_md$partition_loss_mu_prior$partition_loss_mu_prior,
      partition_loss_sigma_prior = pior_md$partition_loss_sigma_prior$partition_loss_sigma_prior,
      partition_loss_max = pior_md$partition_loss_max,
      n_draws = n_draws,
      show_draws = show_draws,
      seed = seed
    )
  }

  return(prior_plot)
}


#' Plot prior distributions for total number of partitions
#'
#' @param max_partitions_prior Prior for the maximum number of partitions.
#' @param partition_loss_mu_prior Prior for the mean of the partition loss.
#' @param partition_loss_sigma_prior Prior for the standard deviation of the
#'   partition loss.
#' @param partition_loss_max Maximum value for the partition loss.
#' @param n_draws Number of draws to simulate from the prior distributions.
#' @param show_draws Number of example draws to show in the plots.
#' @param seed Seed for random number generation to ensure reproducibility.
#'
#' @return A grid of plots showing the distribution of partition loss, 95%
#'   intervals for valid partitions, and distributions of mean and standard
#'   deviation of valid partitions.
#' @keywords internal
plot_prior_partitions <- function(max_partitions_prior, partition_loss_mu_prior,
                                  partition_loss_sigma_prior, partition_loss_max,
                                  n_draws = 1000, show_draws = 50, seed = 0) {
  set.seed(seed)
  if (max_partitions_prior[2] > 0) {
    max_partitions <- extraDistr::rtnorm(
      n_draws, mean = max_partitions_prior[1],
      sd = max_partitions_prior[2], a = 0
    )
  } else {
    max_partitions <- rep(max_partitions_prior[1], n_draws)
  }

  if (partition_loss_mu_prior[2] > 0) {
    partition_loss_mu <- rnorm(
      n_draws, mean = partition_loss_mu_prior[1], sd = partition_loss_mu_prior[2]
    )
  } else {
    partition_loss_mu <- rep(partition_loss_mu_prior[1], n_draws)
  }

  if (partition_loss_sigma_prior[2] > 0) {
    partition_loss_sigma <- extraDistr::rtnorm(
      n_draws, mean = partition_loss_sigma_prior[1],
      sd = partition_loss_sigma_prior[2], a = 0
    )
  } else {
    partition_loss_sigma <- rep(partition_loss_sigma_prior[1], n_draws)
  }

  partition_number_draws <- rbindlist(lapply(1:n_draws, function(i) {
    partition_loss <- partition_loss_max * plogis(
      rnorm(1000, partition_loss_mu[i], partition_loss_sigma[i])
    )
    partition_numbers <- max_partitions[i] * (1 - partition_loss)
    data.table(id = i, loss = partition_loss, partitions = partition_numbers)
  }))

  # compute kernel density estimate for partition number draws across all ids
  max_density <- max(density(partition_number_draws$loss)$y)

  plot_density <- ggplot(partition_number_draws, aes(x=loss)) +
    geom_density(
      data = partition_number_draws[id <= show_draws,],
      linewidth = 0.2, aes(color = as.factor(id))
      ) +
    geom_density(linewidth = 1, color = "black", linetype = "dashed") +
    xlab("Partitions lost (across dPCR runs)") + ylab("Density") +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_bw() +
    coord_cartesian(ylim = c(0, max_density*2)) +
    theme(legend.position = "none") +
    ggtitle("Distribution of partition loss")

  intervals <- partition_number_draws[id <= show_draws, .(
    lower95_partitions = quantile(partitions, probs = 0.025),
    upper95_partitions = quantile(partitions, probs = 0.975)
    ), by = id]
  intervals[, width := upper95_partitions - lower95_partitions]
  setorderv(intervals, "width", order = -1)
  intervals[, order := 1:.N]
  plot_intervals <- ggplot(intervals) +
    geom_errorbar(aes(xmin=lower95_partitions, xmax=upper95_partitions, y=order), linewidth = 0.2, color = "black") +
    xlab("Valid partitions (across dPCR runs)") + ylab("Draw (ordered by width)") +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle("95% intervals for valid partitions")

  # plot distribution of mean partition number across ids
  plot_means <- partition_number_draws[, .(mean_partitions = mean(partitions)), by = id] |>
    ggplot(aes(x=mean_partitions)) +
      geom_histogram(bins = 30, fill = "darkblue", color = "black", alpha = 0.7) +
      xlab("Mean valid partitions") + ylab("Count") +
      theme_bw() +
      theme(legend.position = "none") +
      ggtitle("Distribution of mean valid partitions")

  # plot distribution of sd of partition number across ids
  plot_sds <- partition_number_draws[, .(sd_partitions = sd(partitions)), by = id] |>
    ggplot(aes(x=sd_partitions)) +
      geom_histogram(bins = 30, fill = "darkblue", color = "black", alpha = 0.7) +
      xlab("Standard deviation of valid partitions") + ylab("Count") +
      theme_bw() +
      theme(legend.position = "none") +
      ggtitle("Distribution of sd of valid partitions")

  return(cowplot::plot_grid(
    plot_density, plot_intervals, plot_means, plot_sds,
    ncol = 2, nrow = 2, labels = ""
  ))
}
