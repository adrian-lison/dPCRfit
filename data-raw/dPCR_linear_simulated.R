library(data.table)

# simulation parameters
n_replicates <- 2 # number of technical replicates
n_measured <- 2 # number of biological replicates

alpha_true <- 0.2
beta_true <- 0.6

set.seed(0)

biomass <- sample(seq(0, 28, by = 0.4), size = 50, replace = FALSE)

# simulate data
set.seed(123)

meta_df <- data.table(
  sample_id = paste0("s_",formatC(1:length(biomass), width = 3, flag = "0")),
  biomass = biomass,
  lambda = alpha_true + beta_true * biomass
  )

measurements_df = setNames(data.table(t(sapply(meta_df$lambda, function (lambda) {
  sim_concentrations_dPCR(
    n_samples = n_measured,
    lambda = lambda,
    total_partitions = 25000,
    conversion = 1.73e-5,
    prePCR_cv = 0.1,
    n_replicates = n_replicates,
    distribution = "gamma"
  )
  }
))), as.character(1:n_measured))

# combine
dPCR_linear_simulated <- cbind(meta_df, measurements_df) |>
  melt(
    id.vars = c("sample_id", "biomass", "lambda"),
    variable.name = "replicate_id",
    value.name = "concentration"
    )
dPCR_linear_simulated[, n_technical_reps := n_replicates]

# reorder columns, sample_id and replicate_id first
setcolorder(dPCR_linear_simulated, c(
  "sample_id", "replicate_id", "biomass",
  "lambda", "n_technical_reps", "concentration"
  ))

# oder by sample_id and replicate_id
setorderv(dPCR_linear_simulated, c("sample_id", "replicate_id"))

usethis::use_data(dPCR_linear_simulated, overwrite = TRUE)
