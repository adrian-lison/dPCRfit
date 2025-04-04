functions {
  #include functions/helper_functions.stan
  #include functions/link.stan
  #include functions/dist_normal.stan
  #include functions/dist_lognormal.stan
  #include functions/dist_gamma.stan
  #include functions/hurdle.stan
  #include functions/pcr_noise.stan
}
data {
  // Measurements ----
  int<lower=0> n_samples; // number of different samples
  int<lower=0> n_measured; // number of all measurements
  array[n_measured] int<lower=1, upper=n_samples> measure_to_sample; // index mapping measurements to samples
  vector<lower=0>[n_measured] measured_concentrations; // measured concentrations
  vector<lower=0>[n_measured] n_averaged; // number of averaged technical replicates per measurement (is vector for vectorization)
  vector<lower=0>[n_measured] dPCR_total_partitions; // total number of partitions in dPCR
  int<lower=0, upper=3> obs_dist; // Parametric distribution for observation likelihood: 0 (default) for gamma, 1 for log-normal, 2 for truncated normal, 3 for normal

  // Concentration ----
  array[2] real alpha_prior; // prior for intercept of concentration
  int<lower=0> K; // number of concentration predictors
  matrix[n_samples, K] X; // predictor design matrix
  array[K > 0 ? 2 : 0] real beta_prior; // prior for regression coefficients
  int <lower=0, upper=1> link_type; // 0 for identity, 1 for log-link

  // Coefficient of variation (CV) of measurements ----
  int<lower=0, upper=2> cv_type; // 0 for constant, 1 for dPCR, 2 for constant_var
  array[2] real nu_upsilon_a_prior; // prior for pre-PCR CV
  int<lower=0, upper=1> total_partitions_observe; // 0 for not observed, 1 for observed
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_b_mu_prior; // prior for parameter 2 of CV formula (number of partitions). Scaled by 1e-4 for numerical efficiency.
  array[cv_type == 1 && total_partitions_observe!=1 ? 2 : 0] real nu_upsilon_b_cv_prior; // prior on partition number variation (coefficient of variation)
  array[cv_type == 1 ? 2 : 0] real nu_upsilon_c_prior; // prior for parameter 3 of CV formula (partition size*(scaling factor, i.e. exp_conc_assay/exp_conc_ww)). Scaled by 1e+5 for numerical efficiency.
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_type; // 0 for gamma, 1 for log-normal
  array[cv_type == 1 ? 1 : 0] int <lower=0, upper=1> cv_pre_approx_taylor; // 0 for no Taylor expansion approximation, 1 for Taylor expansion approximation

  // Limit of detection ----
  // LOD_model = 0: no LOD
  // LOD_model = 1: assumed LOD, LOD_scale provided
  // LOD_model = 2: estimated LOD based on dPCR model, needs dPCR parameters
  int<lower=0, upper=2> LOD_model;
  array[(LOD_model == 1) ? 1 : 0] real<lower=0> LOD_scale;
  real<lower=0, upper=1> LOD_drop_prob; // probability threshold for non-detection below which log likelihood contributions of observed concentrations are dropped from LOD model
}
transformed data {

  // number of averaged technical replicates per date
  real n_averaged_median = quantile(n_averaged, 0.5);
  vector[n_samples] n_averaged_all = rep_vector(n_averaged_median, n_samples);
  for (i in 1:n_measured) {
    // note that if several measurements per sample exist,
    // the number of replicates of the last one will be used for that date
    n_averaged_all[measure_to_sample[i]] = n_averaged[i];
  }

  // number of total partitions per measurement per date
  real total_partitions_median = quantile(dPCR_total_partitions, 0.5);
  vector[n_samples] total_partitions_all = rep_vector(total_partitions_median, n_samples);
  for (i in 1:n_measured) {
    // note that if several measurements per sample exist,
    // the total partitions of the last one will be used for that date
    total_partitions_all[measure_to_sample[i]] = dPCR_total_partitions[i];
  }

  // Upper relevant bound for LOD model
  real conc_drop_prob;
  if (LOD_model == 0) {
    conc_drop_prob = positive_infinity();
  } else {
    real LOD_expected_scale;
    if (LOD_model == 1) {
      LOD_expected_scale = LOD_scale[1];
    } else if (LOD_model == 2) {
      LOD_expected_scale = (
      (total_partitions_observe ? total_partitions_median : (nu_upsilon_b_mu_prior[1] * 1e4)) *
      nu_upsilon_c_prior[1] * 1e-5 *
      n_averaged_median
      );
    }
    conc_drop_prob = -log(LOD_drop_prob)/LOD_expected_scale; // concentrations above this value are irrelevant for LOD model (probability of non-detection is virtually zero)
  }

  int n_zero = num_zeros(measured_concentrations);
  int n_dropLOD = num_zeros(fmax(0, conc_drop_prob - measured_concentrations));
  array[n_zero] int<lower=0> i_zero;
  array[n_measured - n_zero] int<lower=0> i_nonzero;
  array[n_measured - n_dropLOD] int<lower=0> i_LOD;
  array[n_measured - n_zero - n_dropLOD] int<lower=0> i_nonzero_LOD;
  int i_z = 0;
  int i_nz = 0;
  int i_lod = 0;
  int i_nzs = 0;
  for (n in 1:n_measured) {
    if (measured_concentrations[n] == 0) {
      i_z += 1;
      i_zero[i_z] = n;
    } else {
      i_nz += 1;
      i_nonzero[i_nz] = n;
      if (measured_concentrations[n] < conc_drop_prob) {
        i_nzs += 1;
        i_nonzero_LOD[i_nzs] = n;
      }
    }
    if (measured_concentrations[n] < conc_drop_prob) {
        i_lod += 1;
        i_LOD[i_lod] = n;
    }
  }
}
parameters {
  // Concentration regression
  real<lower=(link_type == 0 ? 0 : negative_infinity())> alpha;
  vector[K] beta;

  // Coefficient of variation of likelihood for measurements
  real<lower=0> nu_upsilon_a;
  array[(cv_type == 1) && (nu_upsilon_b_mu_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_mu;
  array[(cv_type == 1) && total_partitions_observe!=1 && (nu_upsilon_b_cv_prior[2] > 0) ? 1 : 0] real<lower=0> nu_upsilon_b_cv;
  vector[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b_noise_raw;
  array[(cv_type == 1) && nu_upsilon_c_prior[2] > 0 ? 1 : 0] real<lower=0> nu_upsilon_c;
}
transformed parameters {
  vector<lower=0>[n_samples] true_concentration;
  vector<lower=0>[(cv_type == 1) && total_partitions_observe!=1 ? n_measured : 0] nu_upsilon_b; // total partitions per measurement
  array[LOD_model > 0 ? 1 : 0] vector<lower=0>[n_measured] LOD_hurdle_scale;
  vector[n_measured] concentration;
  vector[n_measured] p_zero_log = rep_vector(negative_infinity(), n_measured);
  vector[n_measured] p_zero = rep_vector(0, n_measured);

  true_concentration = alpha + X * beta;
  if (link_type == 1) {
    true_concentration = exp(true_concentration);
  }
  concentration = true_concentration[measure_to_sample];

  if ((cv_type == 1) && total_partitions_observe!=1) {
    nu_upsilon_b = total_partitions_noncentered(
      param_or_fixed(nu_upsilon_b_mu, nu_upsilon_b_mu_prior),
      param_or_fixed(nu_upsilon_b_cv, nu_upsilon_b_cv_prior),
      nu_upsilon_b_noise_raw
    );
  }

 if (LOD_model > 0) {
    if (LOD_model == 1) {
      LOD_hurdle_scale[1] = rep_vector(LOD_scale[1], n_measured);
    } else if (LOD_model == 2) {
      LOD_hurdle_scale[1] = (
      n_averaged .*
      (total_partitions_observe ? dPCR_total_partitions : nu_upsilon_b * 1e4) *
       param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5
      );
    }
    p_zero_log[i_LOD] = log_hurdle_exponential(
        concentration[i_LOD],
        LOD_hurdle_scale[1][i_LOD], // LOD scale (c * m * n)
        cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
        cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
        );
    p_zero = exp(p_zero_log);
  }

    // CV of each observation as a function of concentration
  vector[n_measured] cv;
  if (cv_type == 0) { // constant cv
    cv = rep_vector(nu_upsilon_a, n_measured);
  } else if (cv_type == 1) { // dPCR
    cv = cv_dPCR_pre(
      concentration, // lambda (concentration)
      nu_upsilon_a, // nu_pre (pre-PCR CV)
      (total_partitions_observe ? dPCR_total_partitions : nu_upsilon_b * 1e4), // m (number of partitions)
      param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5, // c (conversion factor)
      n_averaged, // n (number of averaged replicates)
      cv_pre_type[1], // Type of pre-PCR CV
      cv_pre_approx_taylor[1] // Should taylor approximation be used?
      );
  } else if (cv_type == 2) { // constant variance
    cv = (
      nu_upsilon_a * mean(measured_concentrations) /
      concentration
      );
  }

  vector[n_measured - n_zero] mean_conditional = concentration[i_nonzero] ./ (1-p_zero[i_nonzero]);
  vector[n_measured - n_zero] cv_conditional = sqrt(trim_or_reject_lb(
      cv[i_nonzero]^2 .* (1-p_zero[i_nonzero]) - p_zero[i_nonzero],
      1e-5, // trim to almost zero
      -1 // throw error when significantly below zero
    ));
  if (is_nan(sum(cv_conditional))) {
    for (i in 1:num_elements(cv_conditional)) {
      if (is_nan(cv_conditional[i])) {
        print("i:", i);
        print("p_zero: ", p_zero[i_nonzero][i]);
        print("concentration: ", concentration[i_nonzero][i]);
        print("CV: ", cv[i_nonzero][i]);
        print("alpha", alpha);
        print("beta", beta);
        print("nu_upsilon_a: ", nu_upsilon_a);
        print("nu_upsilon_b: ", nu_upsilon_b * 1e4);
        print("nu_upsilon_c: ", param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5);
      }
    }
  }
}
model {
  // Priors

  // True concentration
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  if (K > 0) {
    beta ~ normal(beta_prior[1], beta_prior[2]);
  }

  // Prior on cv of likelihood for measurements
  nu_upsilon_a ~ normal(nu_upsilon_a_prior[1], nu_upsilon_a_prior[2]) T[0, ]; // truncated normal
  if (cv_type == 1) {
    target += normal_prior_lpdf(nu_upsilon_b_mu | nu_upsilon_b_mu_prior, 0); // truncated normal
    target += normal_prior_lpdf(nu_upsilon_c | nu_upsilon_c_prior, 0); // truncated normal
    if (total_partitions_observe != 1) {
      target += normal_prior_lpdf(nu_upsilon_b_cv | nu_upsilon_b_cv_prior, 0); // truncated normal
      nu_upsilon_b_noise_raw ~ std_normal();
    }
  }

  // Likelihood
  {
    // limit of detection
    if (LOD_model > 0) {
      target += sum(p_zero_log[i_zero]); // below-LOD probabilities
      target += sum(log1m_exp(p_zero_log[i_nonzero_LOD])); // above-LOD probabilities
    }

    // measurements
    if (obs_dist == 0) {
      target += gamma3_lpdf(
        measured_concentrations[i_nonzero] |
        mean_conditional, // expectation
        cv_conditional // coefficient of variation
      );
    } else if (obs_dist == 1) {
      target += lognormal5_lpdf(
        measured_concentrations[i_nonzero] |
        mean_conditional, // log expectation
        cv_conditional // coefficient of variation
      );
    } else if (obs_dist == 2) {
      target += normal2_lpdf(
        measured_concentrations[i_nonzero] |
        mean_conditional, // expectation
        cv_conditional, // coefficient of variation,
        0 // truncate at zero
      );
    } else if (obs_dist == 3) {
      target += normal2_lpdf(
        measured_concentrations |
        concentration, // expectation
        cv // coefficient of variation
      );
    } else {
      reject("Distribution not supported.");
    }
  }
}
generated quantities {
  vector[n_measured] predicted_concentration;
  vector[n_measured] residuals;

  {
    vector[n_samples] nu_upsilon_b_all;
    if (cv_type == 1) {
      if (total_partitions_observe == 1) {
        nu_upsilon_b_all = total_partitions_all / 1e4;
      } else {
        nu_upsilon_b_all = total_partitions_noncentered(
          param_or_fixed(nu_upsilon_b_mu, nu_upsilon_b_mu_prior),
          param_or_fixed(nu_upsilon_b_cv, nu_upsilon_b_cv_prior),
          std_normal_n_rng(n_samples)
          );
        for (i in 1:n_measured) {
          nu_upsilon_b_all[measure_to_sample[i]] = nu_upsilon_b[i];
        }
      }
    }

    vector[n_samples] p_zero_all;
    if (LOD_model > 0) {
      vector[n_samples] LOD_hurdle_scale_all;
      // scale for LOD hurdle model
      if (LOD_model == 1) {
        LOD_hurdle_scale_all = rep_vector(LOD_scale[1], n_samples);
      } else if (LOD_model == 2) {
        LOD_hurdle_scale_all = (
        n_averaged_all .*
        nu_upsilon_b_all * 1e4 *
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5
        );
      }
      p_zero_all = exp(log_hurdle_exponential(
          true_concentration, // lambda (concentration)
          LOD_hurdle_scale_all,
          cv_type == 1 ? nu_upsilon_a : 0, // nu_pre (pre-PCR CV)
          cv_type == 1 ? cv_pre_type[1] : 0 // Type of pre-PCR CV
          ));
      p_zero_all = trim_or_reject_ub(
        p_zero_all,
        1-1e-5, // trim to almost 1
        1.01 // throw error when significantly above 1
      );
    } else {
      p_zero_all = rep_vector(0, n_samples);
    }

    vector[n_samples] cv_all;
    if (cv_type == 0) {
      cv_all = rep_vector(nu_upsilon_a, n_samples);
    } else if (cv_type == 1) {
      cv_all = cv_dPCR_pre(
        true_concentration, // lambda (concentration)
        nu_upsilon_a, // nu_pre (pre-PCR CV)
        nu_upsilon_b_all * 1e4, // m (number of partitions)
        param_or_fixed(nu_upsilon_c, nu_upsilon_c_prior) * 1e-5, // c (conversion factor)
        n_averaged_all, // n (number of averaged replicates) for all dates
        cv_pre_type[1], // Type of pre-PCR CV
        cv_pre_approx_taylor[1] // Should taylor approximation be used?
        );
    } else if (cv_type == 2) {
      cv_all = (
        nu_upsilon_a * mean(measured_concentrations[i_nonzero]) /
        true_concentration
        );
    }

    vector[n_measured] mean_conditional_all = (true_concentration ./ (1-p_zero_all))[measure_to_sample];
    vector[n_measured] cv_conditional_all = sqrt(trim_or_reject_lb(
      cv_all^2 .* (1-p_zero_all) - p_zero_all,
      1e-5, // trim to almost zero
      -1 // throw error when significantly below zero
    ))[measure_to_sample];

    vector[n_measured] isnonzero;
    if (LOD_model > 0) {
      isnonzero = to_vector(bernoulli_rng(1-p_zero_all[measure_to_sample]));
    } else {
      isnonzero = rep_vector(1, n_measured);
    }

    vector[n_measured] meas_conc;
    if (obs_dist == 0) {
      meas_conc = gamma3_rng(mean_conditional_all, cv_conditional_all);
    } else if (obs_dist == 1) {
      meas_conc = lognormal5_rng(mean_conditional_all, cv_conditional_all);
    } else if (obs_dist == 2) {
      meas_conc = normal2_rng(mean_conditional_all, cv_conditional_all, 0); // truncated at zero
    } else if (obs_dist == 3) {
      meas_conc = normal2_rng(true_concentration[measure_to_sample], cv_conditional_all);
    } else {
      reject("Distribution not supported.");
    }
    predicted_concentration = isnonzero .* meas_conc;
    residuals = predicted_concentration - measured_concentrations;
  }
}
