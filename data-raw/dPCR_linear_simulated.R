library(data.table)

# simulation parameters
n_replicates <- 2 # number of technical replicates
n_measured <- 2 # number of biological replicates

alpha_true <- 0.2
beta_true <- 0.6

set.seed(0)

biomass <- sample(seq(0, 28, by = 0.4), size = 50, replace = FALSE)

# simulate data
set.seed(42)

meta_df <- data.table(
  sample_id = paste0("s_",formatC(1:length(biomass), width = 3, flag = "0")),
  biomass = biomass,
  lambda = alpha_true + beta_true * biomass
  )

dPCR_linear_simulated = rbindlist(apply(meta_df, 1, function (meta_df_row) {
  df <- sim_concentrations_dPCR(
    n_samples = n_measured,
    lambda = as.numeric(meta_df_row[["lambda"]]),
    m_partitions = sim_partitions(
      n_measured, n_replicates,
      max_partitions = 30000,
      partition_loss_logit_mean = -2,
      partition_loss_logit_sd = 0.4,
      partition_loss_max = 0.5
    ),
    conversion = 1.73e-5,
    prePCR_cv = 0.1,
    n_replicates = n_replicates,
    distribution = "gamma",
    get_partitions = TRUE
  )
  setDT(df)
  df[, sample_id := meta_df_row[["sample_id"]]]
  df[, biomass := as.numeric(meta_df_row[["biomass"]])]
  df[, lambda := as.numeric(meta_df_row[["lambda"]])]
  df[, bio_replicate_id := 1:.N]
  df[, n_technical_reps := n_replicates]
  return(df)
  }
))

setcolorder(dPCR_linear_simulated, neworder = c(
  "sample_id", "biomass", "lambda", "bio_replicate_id", "n_technical_reps",
  "total_partitions", "positive_partitions", "concentration"
))

data.table::setnames(
  dPCR_linear_simulated,
  c("total_partitions", "positive_partitions"),
  c("avg_total_partitions", "avg_positive_partitions"),
  )

# oder by sample_id and replicate_id
setorderv(dPCR_linear_simulated, c("sample_id", "bio_replicate_id"))

usethis::use_data(dPCR_linear_simulated, overwrite = TRUE)
