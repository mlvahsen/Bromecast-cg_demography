## Preliminaries ####
# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork); library(lubridate); 
library(loo)
# Source data for modeling -- Becca you can ignore this
# source("supp_code/data_prep.R")

# Read in data
cg_model <- read_csv("supp_data/model_data.csv")

## FECUNDITY MODEL WITH RANDOM SLOPES ####
# Make data set of just plants that survived and reproduced
cg_model$survived <- ifelse((cg_model$inflor_mass) > 0, 1, 0)

# Make subset of data that survived
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

# Set contrasts
options(contrasts = c("contr.sum", "contr.poly"))

# Center and scale continuous variables, make factors, make survival response
# binary
cg_model_fecun$clim_dist_sc <- scale(cg_model_fecun$clim_dist)[,1]
# Get quadratic of non-absolute climate distance
cg_model_fecun$clim_dist_sc2 <- (cg_model_fecun$clim_dist_sc)^2
cg_model_fecun$sqrt_new_neighbors_sc <- scale(sqrt(cg_model_fecun$new_neighbors))[,1]
cg_model_fecun$temp_fecun_sc <- scale(cg_model_fecun$temp_fecun)[,1]
cg_model_fecun$vwc_avg_sc <- scale(sqrt(cg_model_fecun$vwc_avg))[,1]
cg_model_fecun$temp_surv_sc <- scale(sqrt(cg_model_fecun$temp_surv))[,1]
cg_model_fecun$genotype <- as.factor(cg_model_fecun$genotype)
cg_model_fecun$survived <- ifelse(cg_model_fecun$survived == "Y", 1, 0)
cg_model_fecun$site_year_gravel <- as.factor(paste(cg_model_fecun$site,
                                                     cg_model_fecun$year,
                                                     cg_model_fecun$albedo, sep = "_"))

# Make design matrix for fixed effects
tibble(nb = cg_model_fecun$sqrt_new_neighbors_sc,
       temp = cg_model_fecun$temp_fecun_sc,
       vwc = cg_model_fecun$vwc_avg_sc,
       clim_dist = cg_model_fecun$clim_dist_sc,
       clim_dist2 = cg_model_fecun$clim_dist_sc2) %>% 
  mutate(nb_temp = nb * temp,
         nb_vw = nb * vwc,
         temp_vwc = temp * vwc) %>% 
  as.matrix() -> X

# Create numeric identifiers for random effects (numeric so the STAN model can
# loop through them; they are still treated as factors)
site_year_id <- as.numeric(cg_model_fecun$site_year_gravel)
plot_unique_id <- as.numeric(as.factor(paste(cg_model_fecun$plot_unique, cg_model_fecun$year)))
genotype_id <- as.numeric(as.factor(cg_model_fecun$genotype))

# Create list of all data needed for the model
data <- list(X = X,
             N = nrow(X),
             P = ncol(X),
             K_site = length(unique(site_year_id)),
             K_plot = length(unique(plot_unique_id)),
             K_genotype = length(unique(genotype_id)),
             site_id = site_year_id,
             plot_id = plot_unique_id,
             genotype_id = genotype_id,
             seed_count = cg_model_fecun$seed_count)

# Settings for STAN run
file <- file.path("supp_code/cmdstanr_write_stan_file_dir/demo_model_fecun.stan")
mod <- cmdstan_model(file, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))
n_cores = parallel::detectCores()
chain = 3
threads_per_chain = ceiling(n_cores/chain)

# Fit model in STAN
fit <- mod$sample(
  data = data,
  chains = chain,
  seed = 4685,
  parallel_chains = 3,
  show_messages = T,
  refresh = 10,
  iter_warmup = 500,
  threads_per_chain = threads_per_chain,
  iter_sampling = 500
)

# Get summary of all parameters
summary = fit$summary()

# Get all posterior draws for parameters
posterior <- fit$draws()


## FECUNDITY MODEL WITHOUT RANDOM SLOPES ####

# Fit the same model without random slopes 
file_noslopes <- file.path("supp_code/cmdstanr_write_stan_file_dir/demo_model_fecun_noslopes.stan")
mod_noslopes <- cmdstan_model(file_noslopes, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))
fit_noslopes <- mod_noslopes$sample(
  data = data,
  chains = chain,
  seed = 4685,
  parallel_chains = 3,
  show_messages = T,
  refresh = 10,
  iter_warmup = 1000,
  threads_per_chain = threads_per_chain,
  iter_sampling = 1000
)

# To save output file
# fit_noslopes$save_output_files("outputs/")

# Compute the LOOs for both models to compare
full_model <- fit$loo(variables = "log_likelihood_values")
reduced_model <- fit_noslopes$loo(variables = "log_likelihood_values")

loo_compare(full_model, reduced_model)

# Get summary of all parameters
summary = fit$summary()

## SURVIVAL MODEL WITH RANDOM SLOPES ####

# Make design matrix for fixed effects
tibble(dens = ifelse(cg_model$density == "hi", -1, 1),
       temp = scale(cg_model$temp_surv)[,1],
       vwc = scale(cg_model$vwc_avg)[,1],
       clim_dist = scale(cg_model$clim_dist)[,1],
       clim_dist2 = scale(cg_model$clim_dist)[,1]^2) %>% 
  mutate(nb_temp = dens * temp,
         nb_vw = dens * vwc,
         temp_vwc = temp * vwc) %>% 
  as.matrix() -> Xs

# Create identifiers for random effects
site_year_ids <- as.numeric(factor(paste(cg_model$site, cg_model$year, cg_model$albedo)))
plot_unique_ids <- as.numeric(as.factor(paste(cg_model$plot_unique, cg_model$year)))
genotype_ids <- as.numeric(as.factor(cg_model$genotype))

# Create data list for survival model
data <- list(X = Xs,
             N = nrow(Xs),
             P = ncol(Xs),
             K_site = length(unique(site_year_ids)),
             K_plot = length(unique(plot_unique_ids)),
             K_genotype = length(unique(genotype_ids)),
             site_id = site_year_ids,
             plot_id = plot_unique_ids,
             genotype_id = genotype_ids,
             survival = cg_model$survived)

# Settings for STAN run
file <- file.path("supp_code/cmdstanr_write_stan_file_dir/demo_model_surv.stan")
mod_s <- cmdstan_model(file, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))

# Fit model in STAN
fit_s <- mod_s$sample(
  data = data,
  chains = chain,
  seed = 4685,
  parallel_chains = 3,
  show_messages = T,
  refresh = 10,
  iter_warmup = 500,
  threads_per_chain = threads_per_chain,
  iter_sampling = 500
)

# To save output files
# fit_s$save_output_files("outputs/")

# To reload output files so we don't have to run model each time
# files = c("outputs/demo_model_surv-202502111510-1-84e4ce.csv",
#           "outputs/demo_model_surv-202502111510-2-84e4ce.csv",
#           "outputs/demo_model_surv-202502111510-3-84e4ce.csv")
# fit_s <- as_cmdstan_fit(files)

## SRUVIVAL MODEL WITHOUT RANDOM SLOPES ####

file_s <- file.path("supp_code/cmdstanr_write_stan_file_dir/demo_model_surv_noslopes.stan")
mod_s_noslopes <- cmdstan_model(file_s, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))

fit_s_noslopes <- mod_s_noslopes$sample(
  data = data,
  chains = chain,
  seed = 4685,
  parallel_chains = 3,
  show_messages = T,
  refresh = 10,
  iter_warmup = 500,
  threads_per_chain = threads_per_chain,
  iter_sampling = 500
)

# Compute the LOOs for both models to compare
full_model_s <- fit_s$loo(variables = "log_likelihood_values")
reduced_model_s <- fit_s_noslopes$loo(variables = "log_likelihood_values")
loo_compare(full_model_s, reduced_model_s)

# Get summary of all parameters
summary_s = fit_s$summary()

# Get all posterior draws for parameters
posterior_s <- fit_s$draws()

