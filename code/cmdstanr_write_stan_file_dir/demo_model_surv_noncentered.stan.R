data {
  int<lower=1> K_site;
  int<lower=1> K_plot;
  int<lower=1> K_genotype;
  int<lower=0> N;
  int<lower=0> P;
  array[N] int<lower=1, upper=K_site> site_id;
  array[N] int<lower=1, upper=K_plot> plot_id;
  array[N] int<lower=1, upper=K_genotype> genotype_id;
  matrix[N, P] X;
  array[N] int<lower=0, upper=1> survival;
}

parameters {
  real alpha;
  vector[P] beta;
  
  real<lower=0> sigma_site;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_genotype;
  
  real<lower=0> sigma_genotype_slope_neighbors;
  real<lower=0> sigma_genotype_slope_temp;
  real<lower=0> sigma_genotype_slope_vwc;
  
  // Standard normal latent variables for non-centered parameterization
  vector[K_site] z_site;
  vector[K_plot] z_plot;
  vector[K_genotype] z_genotype;
  matrix[K_genotype, 3] z_genotype_slopes;
}

transformed parameters {
  // Non-centered parameterization
  vector[K_site] site_intercepts = z_site * sigma_site;
  vector[K_plot] plot_intercepts = z_plot * sigma_plot;
  vector[K_genotype] genotype_intercepts = z_genotype * sigma_genotype;
  
  matrix[K_genotype, 3] genotype_slopes;
  for (g in 1:K_genotype) {
    genotype_slopes[g, 1] = z_genotype_slopes[g, 1] * sigma_genotype_slope_neighbors;
    genotype_slopes[g, 2] = z_genotype_slopes[g, 2] * sigma_genotype_slope_temp;
    genotype_slopes[g, 3] = z_genotype_slopes[g, 3] * sigma_genotype_slope_vwc;
  }
}

model {
  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  
  sigma_site ~ normal(0, 1);
  sigma_plot ~ normal(0, 1);
  sigma_genotype ~ normal(0, 1);
  
  sigma_genotype_slope_neighbors ~ normal(0, 1);
  sigma_genotype_slope_temp ~ normal(0, 1);
  sigma_genotype_slope_vwc ~ normal(0, 1);
  
  // Standard normal priors for latent variables
  z_site ~ normal(0, 1);
  z_plot ~ normal(0, 1);
  z_genotype ~ normal(0, 1);
  to_vector(z_genotype_slopes) ~ normal(0, 1);
  
  // Likelihood
  for (n in 1:N) {
    real logit_p = alpha + dot_product(X[n], beta)
    + site_intercepts[site_id[n]]
    + plot_intercepts[plot_id[n]]
    + genotype_intercepts[genotype_id[n]]
    + genotype_slopes[genotype_id[n], 1] * X[n, 1]  // neighbors
    + genotype_slopes[genotype_id[n], 2] * X[n, 2]  // temp
    + genotype_slopes[genotype_id[n], 3] * X[n, 3]; // vwc
    
    survival[n] ~ bernoulli_logit(logit_p);
  }
}

generated quantities {
  vector[N] survival_pred;
  vector[N] log_likelihood_values;
  
  for (n in 1:N) {
    real logit_p = alpha + dot_product(X[n], beta)
    + site_intercepts[site_id[n]]
    + plot_intercepts[plot_id[n]]
    + genotype_intercepts[genotype_id[n]]
    + genotype_slopes[genotype_id[n], 1] * X[n, 1]
    + genotype_slopes[genotype_id[n], 2] * X[n, 2]
    + genotype_slopes[genotype_id[n], 3] * X[n, 3];
    
    survival_pred[n] = bernoulli_logit_rng(logit_p);
    log_likelihood_values[n] = bernoulli_logit_lpmf(survival[n] | logit_p);
  }
}
