## Preliminaries ####
# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork); library(lubridate); 
library(loo); library(cmdstanr)

# Source data for modeling 
#source("supp_code/data_prep.R")
cg_model <- read_csv("data/cg_model_data.csv")

# Make data set of just plants that survived and reproduced
cg_model$survived <- ifelse(cg_model$seed_count > 0, 1, 0)

# Make subset of data that survived
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

# Center and scale continuous variables, make factors, make survival response
# binary
cg_model_fecun$clim_dist_sc <- scale(cg_model_fecun$clim_dist)[,1]
# Get quadratic of climate difference
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
## FECUNDITY MODEL WITH RANDOM SLOPES ####

# Make design matrix for fixed effects
tibble(nb = cg_model_fecun$sqrt_new_neighbors_sc,
       temp = cg_model_fecun$temp_fecun_sc,
       vwc = cg_model_fecun$vwc_avg_sc,
       clim_dist = cg_model_fecun$clim_dist_sc, 
       clim_dist2 = cg_model_fecun$clim_dist_sc2) %>% 
  mutate(nb_temp = nb * temp,
         nb_vwc = nb * vwc,
         climdist_nb = clim_dist * nb,
         climdist2_nb = clim_dist2 * nb) %>% 
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
  iter_warmup = 1000,
  threads_per_chain = threads_per_chain,
  iter_sampling = 2000
)

# Get summary of all parameters
summary = fit$summary()

# Get all posterior draws for parameters
posterior <- fit$draws(format = "df")

# To save output files
fit$save_output_files("outputs/")
  

## FECUNDITY MODEL WITHOUT RANDOM SLOPES ####

# Make design matrix for fixed effects -- need to drop out clim dist information
# too
tibble(nb = cg_model_fecun$sqrt_new_neighbors_sc,
       temp = cg_model_fecun$temp_fecun_sc,
       vwc = cg_model_fecun$vwc_avg_sc) %>% 
  mutate(nb_temp = nb * temp,
         nb_vwc = nb * vwc) %>% 
  as.matrix() -> X_noclimdist

# Create list of all data needed for the model
data_noslopes <- list(X = X_noclimdist,
             N = nrow(X_noclimdist),
             P = ncol(X_noclimdist),
             K_site = length(unique(site_year_id)),
             K_plot = length(unique(plot_unique_id)),
             K_genotype = length(unique(genotype_id)),
             site_id = site_year_id,
             plot_id = plot_unique_id,
             genotype_id = genotype_id,
             seed_count = cg_model_fecun$seed_count)

# Fit the same model without random slopes 
file_noslopes <- file.path("supp_code/cmdstanr_write_stan_file_dir/demo_model_fecun_noslopes.stan")
mod_noslopes <- cmdstan_model(file_noslopes, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))
fit_noslopes <- mod_noslopes$sample(
  data = data_noslopes,
  chains = chain,
  seed = 4685,
  parallel_chains = 3,
  show_messages = T,
  refresh = 10,
  iter_warmup = 1000,
  threads_per_chain = threads_per_chain,
  iter_sampling = 2000
)

# To save output file
fit_noslopes$save_output_files("outputs/")

## Compare fecundity models ####

# Read in models if not running them above
fit <- as_cmdstan_fit(files = c("outputs/demo_model_fecun-202506031125-1-391f4a.csv",
                                "outputs/demo_model_fecun-202506031125-2-391f4a.csv",
                                "outputs/demo_model_fecun-202506031125-3-391f4a.csv"))

fit_noslopes <- as_cmdstan_fit(files = c("outputs/demo_model_fecun_noslopes-202506041129-1-44ac8c.csv",
                                         "outputs/demo_model_fecun_noslopes-202506041129-2-44ac8c.csv",
                                         "outputs/demo_model_fecun_noslopes-202506041129-3-44ac8c.csv"))

# Compute the LOOs for both models to compare
full_model <- fit$loo(variables = "log_likelihood_values")
reduced_model <- fit_noslopes$loo(variables = "log_likelihood_values")

# Info for Table 1
loo_compare(full_model, reduced_model)

fit$draws(format = "df") %>% 
  dplyr::select(starts_with("sigma")) %>% 
  gather(key = parameter, value = value, sigma_genotype:sigma_genotype_slope_vwc) %>% 
  group_by(parameter) %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>% 
  arrange(-mean)

# Get summary of all parameters
summary = fit$summary()

## SURVIVAL MODEL WITH RANDOM SLOPES ####

# Make design matrix for fixed effects
tibble(dens = ifelse(cg_model$density == "hi", -1, 1),
       temp = scale(cg_model$temp_surv)[,1],
       vwc = scale(cg_model$vwc_avg)[,1],
       clim_dist = scale(cg_model$clim_dist)[,1],
       clim_dist2 = scale(cg_model$clim_dist)[,1]^2) %>% 
  mutate(dens_temp = dens * temp,
         dens_vwc = dens * vwc,
         climdist_dens = clim_dist * dens,
         climdist2_dens = clim_dist2 * dens) %>% 
  as.matrix() -> Xs

# Create identifiers for random effects
site_year_ids <- as.numeric(factor(paste(cg_model$site, cg_model$year, cg_model$albedo)))
plot_unique_ids <- as.numeric(as.factor(paste(cg_model$plot_unique, cg_model$year)))
genotype_ids <- as.numeric(as.factor(cg_model$genotype))

# Create data list for survival model
data_s <- list(X = Xs,
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
file <- file.path("supp_code/cmdstanr_write_stan_file_dir/demo_model_surv_noncenter.stan")
mod_s <- cmdstan_model(file, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))

# Fit model in STAN
fit_s <- mod_s$sample(
  data = data_s,
  chains = chain,
  seed = 4685,
  parallel_chains = 3,
  show_messages = T,
  refresh = 10,
  iter_warmup = 1000,
  threads_per_chain = threads_per_chain,
  iter_sampling = 2000
)

# To save output files
fit_s$save_output_files("outputs/")

## SURVIVAL MODEL WITHOUT RANDOM SLOPES ####

# Make design matrix for fixed effects -- drop clim dist as well
tibble(dens = ifelse(cg_model$density == "hi", -1, 1),
       temp = scale(cg_model$temp_surv)[,1],
       vwc = scale(cg_model$vwc_avg)[,1]) %>% 
  mutate(dens_temp = dens * temp,
         dens_vwc = dens * vwc) %>% 
  as.matrix() -> Xs_noclimdist

# Create data list for survival model
data_s_noslopes <- list(X = Xs_noclimdist,
                        N = nrow(Xs_noclimdist),
                        P = ncol(Xs_noclimdist),
                        K_site = length(unique(site_year_ids)),
                        K_plot = length(unique(plot_unique_ids)),
                        K_genotype = length(unique(genotype_ids)),
                        site_id = site_year_ids,
                        plot_id = plot_unique_ids,
                        genotype_id = genotype_ids,
                        survival = cg_model$survived)

file_s <- file.path("supp_code/cmdstanr_write_stan_file_dir/demo_model_surv_noslopes.stan")
mod_s_noslopes <- cmdstan_model(file_s, stanc_options = list("O1"), cpp_options = list(stan_threads = TRUE))

fit_s_noslopes <- mod_s_noslopes$sample(
  data = data_s_noslopes,
  chains = chain,
  seed = 4685,
  parallel_chains = 3,
  show_messages = T,
  refresh = 10,
  iter_warmup = 1000,
  threads_per_chain = threads_per_chain,
  iter_sampling = 2000
)

fit_s_noslopes$save_output_files("outputs/")

# Info for Table 1

# Compute the LOOs for both models to compare
full_model_s <- fit_s$loo(variables = "log_likelihood_values")
reduced_model_s <- fit_s_noslopes$loo(variables = "log_likelihood_values")
loo_compare(full_model_s, reduced_model_s)

fit_s$draws(format = "df") %>% 
  dplyr::select(starts_with("sigma")) %>% 
  gather(key = parameter, value = value, sigma_genotype:sigma_genotype_slope_vwc) %>% 
  group_by(parameter) %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>% 
  arrange(-mean)

# Get summary of all parameters
summary_s = fit_s$summary()

# Get all posterior draws for parameters
posterior_s <- fit_s$draws()

# Read in models if not running them above
fit_s <- as_cmdstan_fit(files = c("outputs/demo_model_surv_noncenter-202506031158-1-25d6a2.csv",
                                "outputs/demo_model_surv_noncenter-202506031158-2-25d6a2.csv",
                                "outputs/demo_model_surv_noncenter-202506031158-3-25d6a2.csv"))

fit_s_noslopes <- as_cmdstan_fit(files = c("outputs/demo_model_surv_noslopes-202506041259-1-68ac02.csv",
                                         "outputs/demo_model_surv_noslopes-202506041259-2-68ac02.csv",
                                         "outputs/demo_model_surv_noslopes-202506041259-3-68ac02.csv"))

full_model_s <- fit_s$loo(variables = "log_likelihood_values")
reduced_model_s <- fit_s_noslopes$loo(variables = "log_likelihood_values")
loo_compare(full_model_s, reduced_model_s)

fit_s$draws(format = "df") %>% 
  dplyr::select(starts_with("sigma")) %>% 
  gather(key = parameter, value = value, sigma_genotype:sigma_genotype_slope_vwc) %>% 
  group_by(parameter) %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>% 
  arrange(-mean)
