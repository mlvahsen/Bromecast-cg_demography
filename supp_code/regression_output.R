## Preliminaries ####
# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork); library(lubridate); 
library(loo); library(cmdstanr)

# Survival model
files = c("outputs/demo_model_surv_noncenter-202506091738-1-51da52.csv",
          "outputs/demo_model_surv_noncenter-202506091738-1-51da52.csv",
          "outputs/demo_model_surv_noncenter-202506091738-1-51da52.csv")
fit_s <- as_cmdstan_fit(files)
draws_s <- fit_s$draws(format = "df")

# Fecundity model
files_f = c("outputs/demo_model_fecun-202506091642-1-6d200d.csv",
            "outputs/demo_model_fecun-202506091642-2-6d200d.csv",
            "outputs/demo_model_fecun-202506091642-3-6d200d.csv")
fit <- as_cmdstan_fit(files_f)
draws_f <- fit$draws(format = "df")

# Get regression coefficients mean and 95% credible intervals
draws_f %>%
  select(contains("beta[")) %>% 
  colMeans() %>% 
  round(digits = 3)
  
draws_f %>%
  select(contains("beta[")) %>% 
  apply(2, quantile, probs = c(0.025, 0.975)) %>%  
  round(digits = 3)

draws_s %>%
  select(contains("beta[")) %>% 
  colMeans() %>% 
  round(digits = 3)

draws_s %>%
  select(contains("beta[")) %>% 
  apply(2, quantile, probs = c(0.025, 0.975)) %>%  
  round(digits = 3)

# Random effects
draws_f %>%
  select(contains("sigma")) %>% 
  colMeans() %>% 
  round(digits = 3)

draws_f %>%
  select(contains("sigma")) %>% 
  apply(2, quantile, probs = c(0.025, 0.975)) %>%  
  round(digits = 3)

draws_s %>%
  select(contains("sigma")) %>% 
  colMeans() %>% 
  round(digits = 3)

draws_s %>%
  select(contains("sigma")) %>% 
  apply(2, quantile, probs = c(0.025, 0.975)) %>%  
  round(digits = 3)

