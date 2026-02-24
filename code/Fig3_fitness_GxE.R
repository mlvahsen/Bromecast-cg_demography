# This code creates Figure 3: G × temperature and soil moisture reaction norms
# and sensitivities for ln(fitness). It also creates FIgures S5 and S6 which
# show the survival and ln(fecundity) component reaction norms and sensitivies
# for temperature and soil moisture, respectively.

## Preliminaries ####

# Load libraries
library(scales); library(tidyverse); library(cmdstanr); library(patchwork)

# Read in model output
fit <- as_cmdstan_fit(files = c("outputs/demo_model_fecun-202506091642-1-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-2-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-3-6d200d.csv"))

fit_s <- as_cmdstan_fit(files = c("outputs/demo_model_surv_noncenter-202506091738-1-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-2-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-3-51da52.csv"))

# Read in data
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
cg_model_fecun$temp_sc <- scale(cg_model_fecun$temp_fecun)[,1]
cg_model_fecun$vwc_avg_sc <- scale(sqrt(cg_model_fecun$vwc_avg))[,1]
cg_model_fecun$temp_surv_sc <- scale(sqrt(cg_model_fecun$temp_surv))[,1]
cg_model_fecun$genotype <- as.factor(cg_model_fecun$genotype)
cg_model_fecun$survived <- ifelse(cg_model_fecun$survived == "Y", 1, 0)
cg_model_fecun$site_year_gravel <- as.factor(paste(cg_model_fecun$site,
                                                   cg_model_fecun$year,
                                                   cg_model_fecun$albedo, sep = "_"))

# Make vector of temperature scaled
seq_sc_temp <- seq(min(cg_model_fecun$temp_sc),
                   max(cg_model_fecun$temp_sc),
                   length.out = 10)
# Make same length vector of temperature unscaled
seq_temp <- seq(min(cg_model_fecun$temp_fecun),
                max(cg_model_fecun$temp_fecun),
                length.out = 10)

# Make vector of soil moisture scaled
seq_sc_vwc <- seq(min(cg_model_fecun$vwc_avg_sc),
                  max(cg_model_fecun$vwc_avg_sc),
                  length.out = 10)
# Make same length vector of soil moisture unscaled
seq_vwc <- seq(min(cg_model_fecun$vwc_avg),
               max(cg_model_fecun$vwc_avg),
               length.out = 10)


## Fecundity subplots (temp) ####

# Extract draws from model fit -- fecundity
draws <- fit$draws(format = "df")

# Collect intercept and temperature regression coefficient
alpha <- draws$alpha
beta_2 <- draws$`beta[2]`

# Precompute the intercepts and slopes for all genotypes (outside of loops)
intercepts <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()

# We assume that genotype_slopes is a matrix with dimensions: (samples x 3*genotypes)
slopes <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",2]")) %>% as.matrix()

# Initialize the output matrix
genotype_store<- matrix(NA, nrow = length(seq_temp), ncol = ncol(intercepts))

genotype_samples_f <- list()

# Vectorized computation of genotype stores
for (j in 1:ncol(intercepts)) {
  # Extract the relevant intercept and slopes for genotype j
  intercept_j <- intercepts[, j]  # Vector of intercepts for genotype j
  
  # Correct indexing to get the slopes for the 2nd predictor for genotype j
  slope_j <- slopes[, j]  # Slopes for the second predictor (2nd column) for genotype j
  
  # Precompute the entire matrix of `store_temp_gs` for each genotype in a vectorized manner
  store_temp_gs <- matrix(NA, nrow = length(seq_temp), ncol = nrow(intercepts))
  
  for (i in 1:length(seq_temp)) {
    store_temp_gs[i, ] <- alpha + intercept_j + (beta_2 + slope_j) * seq_sc_temp[i]
  }
  
  # Apply colMeans to reduce the result
  genotype_store[, j] <- rowMeans(store_temp_gs)
  genotype_samples_f[[j]] <- store_temp_gs
}

as_tibble(genotype_store) -> genotype_store_df
names(genotype_store_df) <- paste("g", 1:ncol(intercepts), sep = "_")
genotype_store_df %>% 
  mutate(temp = seq_temp) %>% 
  gather(key = genotype, value = ln_seed_count, g_1:g_96) %>% 
  mutate(genotype = as.factor(parse_number(genotype))) -> temp_g_pred

# Get genotype numeric codes that were used in fitting the model and line them
# back up with true codes
genotype_id <- as.numeric(as.factor(cg_model_fecun$genotype))

tibble(genotype = factor(genotype_id),
       genotype_id = cg_model_fecun$genotype) %>% 
  distinct() %>% 
  merge(temp_g_pred) %>% 
  dplyr::select(genotype = genotype_id, temp, ln_seed_count) ->temp_g_pred

# Collect PC1 values for each genotype
cg_model %>% 
  dplyr::select(genotype, PC1) %>% 
  distinct() %>% 
  merge(temp_g_pred) -> temp_g_pred

temp_g_pred %>% 
  ggplot(aes(x = temp, y = ln_seed_count, group = genotype, color = PC1)) +
  geom_line(linewidth = 0.5) +
  theme_classic(base_size = 16) +
  labs(y = "ln(seed count)",
       x = "soil temperature (°C)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1)  -> a_temp

temp_g_pred %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = ln_seed_count[which(temp == min(temp))],
          max = ln_seed_count[which(temp == max(temp))],
          slope = (max - min) / (max(temp) - min(temp))) -> slope_fecun

coef(summary(lm(slope ~ PC1, data = slope_fecun)))[2,1] -> slope
confint(lm(slope ~ PC1, data = slope_fecun))[2,]-> conf_int

# Set size for text annotation
beta_size <- 4

slope_fecun %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "ln(seed count)/", Delta, "temp."))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                     limits = c(-0.16, 0.03)) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = 0), color = "gray47", linewidth = 1.2) +
  ylim(-0.15, 0.02)+
  # Adding annotation with slope and 95% CI
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf(
      "italic(beta)==%.3f*' [%.3f, %.3f]'",
      slope, conf_int[1], conf_int[2]
    ),
    parse = TRUE,
    hjust = 1.1, vjust = 1.1,
    size = beta_size
  )-> b_temp

## Survival subplots (temp) ####
# Extract draws from model fit -- survival
draws_s <- fit_s$draws(format = "df")

cg_model$temp_surv_sc <- scale(cg_model$temp_surv)[,1]

seq_sc_temps <- seq(min(cg_model$temp_surv_sc),
                    max(cg_model$temp_surv_sc),
                    length.out = 10)
seq_temps <- seq(min(cg_model$temp_surv),
                 max(cg_model$temp_surv),
                 length.out = 10)

alpha <- draws_s$alpha
beta_2 <- draws_s$`beta[2]`

# Precompute the intercepts and slopes for all genotypes (outside of loops)
intercepts <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()

# We assume that genotype_slopes is a matrix with dimensions: (samples x 3*genotypes)
slopes <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",2]")) %>% as.matrix()

# Initialize the output matrix
genotype_stores <- matrix(NA, nrow = length(seq_sc_temps), ncol = ncol(intercepts))

genotype_samples_s <- list()

# Vectorized computation of genotype stores
for (j in 1:ncol(intercepts)) {
  # Extract the relevant intercept and slopes for genotype j
  intercept_j <- intercepts[, j]  # Vector of intercepts for genotype j
  
  # Correct indexing to get the slopes for the 2nd predictor for genotype j
  slope_j <- slopes[, j]  # Slopes for the second predictor (2nd column) for genotype j
  
  # Precompute the entire matrix of `store_temp_gs` for each genotype in a vectorized manner
  store_temp_gs <- matrix(NA, nrow = length(seq_sc_temp), ncol = nrow(intercepts))
  
  for (i in 1:length(seq_sc_temp)) {
    store_temp_gs[i, ] <- alpha + intercept_j + (beta_2 + slope_j) * seq_sc_temps[i]
  }
  
  # Apply colMeans to reduce the result
  genotype_stores[, j] <- rowMeans(store_temp_gs)
  genotype_samples_s[[j]] <- store_temp_gs
}

as_tibble(genotype_stores) -> genotype_store_dfs
names(genotype_store_dfs) <- paste("g", 1:ncol(intercepts), sep = "_")
genotype_store_dfs %>% 
  mutate(temp_surv = seq_temps) %>% 
  gather(key = genotype, value = logit_survival, g_1:g_96) %>% 
  mutate(genotype = as.factor(parse_number(genotype)),
         survival = plogis(logit_survival)) -> temp_g_preds

# Get genotype numeric codes that were used in fitting the model and line them
# back up with true codes
genotype_ids <- as.numeric(as.factor(cg_model$genotype))

tibble(genotype = factor(genotype_ids),
       genotype_id = cg_model$genotype) %>% 
  distinct() %>% 
  merge(temp_g_preds) %>% 
  dplyr::select(genotype = genotype_id, temp_surv, survival, logit_survival) ->temp_g_preds

# Get in PC1 values by genotype
cg_model %>% 
  dplyr::select(genotype, PC1) %>% 
  distinct() %>% 
  merge(temp_g_preds) -> temp_g_preds

temp_g_preds %>% 
  ggplot(aes(x = temp_surv, y = survival, group = genotype, color = PC1)) +
  geom_line(linewidth = 0.5) +
  theme_classic(base_size = 16) +
  labs(y = "P(survival)",
       x = "soil temperature (°C)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) -> a_temp_s

temp_g_preds %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = logit_survival[which(temp_surv == min(temp_surv))],
          max = logit_survival[which(temp_surv == max(temp_surv))],
          slope = (max - min) / (max(temp_surv) - min(temp_surv))) ->slope_surv
  
coef(summary(lm(slope ~ PC1, data = slope_surv)))[2,1]-> slopes
confint(lm(slope ~ PC1, data = slope_surv))[2,] -> conf_ints


slope_surv %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "P(survival)/", Delta, "temp."))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +
  #geom_hline(aes(yintercept = 0), color = "gray47", linewidth = 1.2) +
  # Adding annotation with slope and 95% CI
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf(
      "italic(beta)==%.3f*' [%.3f, %.3f]'",
      slopes, conf_ints[1], conf_ints[2]
    ),
    parse = TRUE,
    hjust = 1.1, vjust = 1.1,
    size = beta_size
  ) +
  scale_y_continuous(limits = c(0.14, 0.30))-> b_temp_s

## Generate Figure S5 ####
a_temp_s + b_temp_s +
  a_temp + b_temp +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        legend.key.width=unit(1.5,"cm")) &
  labs(color = "PC 1 (hot & dry → cool & wet)") &
  guides(colour = guide_colourbar(title.position="bottom", title.hjust = 0.5)) -> temp_plot

png("figs/FigS5_temp_rxn.png", height = 8, width = 8, res = 300, units = "in")
temp_plot
dev.off()


## Fecundity subplots (vwc) ####
draws <- fit$draws(format = "df")

# Collect intercept and temperature regression coefficient
alpha <- draws$alpha
beta_3 <- draws$`beta[3]`

# Precompute the intercepts and slopes for all genotypes (outside of loops)
intercepts <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()

# We assume that genotype_slopes is a matrix with dimensions: (samples x 3*genotypes)
slopes <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",3]")) %>% as.matrix()

# Initialize the output matrix
genotype_store<- matrix(NA, nrow = length(seq_vwc), ncol = ncol(intercepts))

genotype_samples_f <- list()

# Vectorized computation of genotype stores
for (j in 1:ncol(intercepts)) {
  # Extract the relevant intercept and slopes for genotype j
  intercept_j <- intercepts[, j]  # Vector of intercepts for genotype j
  
  # Correct indexing to get the slopes for the 2nd predictor for genotype j
  slope_j <- slopes[, j]  # Slopes for the second predictor (2nd column) for genotype j
  
  # Precompute the entire matrix of `store_temp_gs` for each genotype in a vectorized manner
  store_vwc_gs <- matrix(NA, nrow = length(seq_vwc), ncol = nrow(intercepts))
  
  for (i in 1:length(seq_vwc)) {
    store_vwc_gs[i, ] <- alpha + intercept_j + (beta_3 + slope_j) * seq_sc_vwc[i]
  }
  
  # Apply colMeans to reduce the result
  genotype_store[, j] <- rowMeans(store_vwc_gs)
  genotype_samples_f[[j]] <- store_vwc_gs
}

as_tibble(genotype_store) -> genotype_store_df
names(genotype_store_df) <- paste("g", 1:ncol(intercepts), sep = "_")
genotype_store_df %>% 
  mutate(vwc_avg = seq_vwc) %>% 
  gather(key = genotype, value = ln_seed_count, g_1:g_96) %>% 
  mutate(genotype = as.factor(parse_number(genotype))) -> vwc_g_pred

# Get genotype numeric codes that were used in fitting the model and line them
# back up with true codes
genotype_id <- as.numeric(as.factor(cg_model_fecun$genotype))

tibble(genotype = factor(genotype_id),
       genotype_id = cg_model_fecun$genotype) %>% 
  distinct() %>% 
  merge(vwc_g_pred) %>% 
  dplyr::select(genotype = genotype_id, vwc_avg, ln_seed_count) ->vwc_g_pred

# Collect PC1 values for each genotype
cg_model %>% 
  dplyr::select(genotype, PC1) %>% 
  distinct() %>% 
  merge(vwc_g_pred) -> vwc_g_pred

vwc_g_pred %>% 
  ggplot(aes(x = vwc_avg, y = ln_seed_count, group = genotype, color = PC1)) +
  geom_line(linewidth = 0.5) +
  #geom_point(data = cg_model_fecun) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "ln(seed count)",
       x = "soil moisture (vwc)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1)  -> a_vwc

vwc_g_pred %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = ln_seed_count[which(vwc_avg == min(vwc_avg))],
          max = ln_seed_count[which(vwc_avg == max(vwc_avg))],
          slope = (max - min) / (max(vwc_avg) - min(vwc_avg))) -> slope_fecun

slope <- coef(lm(slope ~ PC1, data = slope_fecun))[2]
conf_int <- confint(lm(slope ~ PC1, data = slope_fecun))[2,]

# Set size for text annotation
beta_size <- 4

slope_fecun %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "ln(seed count)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +
  # Adding annotation with slope and 95% CI
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf(
      "italic(beta)==%.3f*' [%.3f, %.3f]'",
      slope, conf_int[1], conf_int[2]
    ),
    parse = TRUE,
    hjust = 1.1, vjust = 1.1,
    size = beta_size
  )-> b_vwc

## Survival subplot (vwc) ####

cg_model$vwc_avg_sc <- scale(cg_model$vwc_avg)[,1]

seq_sc_vwc <- seq(min(cg_model$vwc_avg_sc),
                  max(cg_model$vwc_avg_sc),
                  length.out = 10)
seq_vwc <- seq(min(cg_model$vwc_avg),
               max(cg_model$vwc_avg),
               length.out = 10)

# Extract draws from model fit -- survival
draws_s <- fit_s$draws(format = "df")

alpha <- draws_s$alpha
beta_3 <- draws_s$`beta[3]`

# Precompute the intercepts and slopes for all genotypes (outside of loops)
intercepts <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()

# We assume that genotype_slopes is a matrix with dimensions: (samples x 3*genotypes)
slopes <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",3]")) %>% as.matrix()

# Initialize the output matrix
genotype_stores <- matrix(NA, nrow = length(seq_sc_vwc), ncol = ncol(intercepts))

genotype_samples_s <- list()

# Vectorized computation of genotype stores
for (j in 1:ncol(intercepts)) {
  # Extract the relevant intercept and slopes for genotype j
  intercept_j <- intercepts[, j]  # Vector of intercepts for genotype j
  
  # Correct indexing to get the slopes for the 2nd predictor for genotype j
  slope_j <- slopes[, j]  # Slopes for the second predictor (2nd column) for genotype j
  
  # Precompute the entire matrix of `store_temp_gs` for each genotype in a vectorized manner
  store_vwc_gs <- matrix(NA, nrow = length(seq_sc_vwc), ncol = nrow(intercepts))
  
  for (i in 1:length(seq_sc_vwc)) {
    store_vwc_gs[i, ] <- alpha + intercept_j + (beta_3 + slope_j) * seq_sc_vwc[i]
  }
  
  # Apply colMeans to reduce the result
  genotype_stores[, j] <- rowMeans(store_vwc_gs)
  genotype_samples_s[[j]] <- store_vwc_gs
}

as_tibble(genotype_stores) -> genotype_store_dfs
names(genotype_store_dfs) <- paste("g", 1:ncol(intercepts), sep = "_")
genotype_store_dfs %>% 
  mutate(vwc_avg = seq_vwc) %>% 
  gather(key = genotype, value = logit_survival, g_1:g_96) %>% 
  mutate(genotype = as.factor(parse_number(genotype)),
         survival = plogis(logit_survival)) -> vwc_g_preds

# Get genotype numeric codes that were used in fitting the model and line them
# back up with true codes
genotype_ids <- as.numeric(as.factor(cg_model$genotype))

tibble(genotype = factor(genotype_ids),
       genotype_id = cg_model$genotype) %>% 
  distinct() %>% 
  merge(vwc_g_preds) %>% 
  dplyr::select(genotype = genotype_id, vwc_avg, survival, logit_survival) ->vwc_g_preds

# Get in PC1 values by genotype
cg_model %>% 
  dplyr::select(genotype, PC1) %>% 
  distinct() %>% 
  merge(vwc_g_preds) -> vwc_g_preds

vwc_g_preds %>% 
  ggplot(aes(x = vwc_avg, y = survival, group = genotype, color = PC1)) +
  geom_line(linewidth = 0.5) +
  #geom_point(data = cg_model_fecun) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "P(survival)",
       x = "soil moisture (vwc)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) -> a_vwc_s

vwc_g_preds %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = logit_survival[which(vwc_avg == min(vwc_avg))],
          max = logit_survival[which(vwc_avg == max(vwc_avg))],
          slope = (max - min) / (max(vwc_avg) - min(vwc_avg))) -> slope_survival

coef(lm(slope ~ PC1, data = slope_survival))[2] -> slopes
confint(lm(slope ~ PC1, data = slope_survival))[2,] -> conf_ints

slope_survival %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "P(survival)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +
  # Adding annotation with slope and 95% CI
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf(
      "italic(beta)==%.3f*' [%.3f, %.3f]'",
      slopes, conf_ints[1], conf_ints[2]
    ),
    parse = TRUE,
    hjust = 1.1, vjust = 1.1,
    size = beta_size
  )-> b_vwc_s

## Generate Figure S6 ####
a_vwc_s + b_vwc_s +
  a_vwc + b_vwc +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        legend.key.width=unit(1.5,"cm")) &
  labs(color = "PC 1 (hot & dry → cool & wet)") &
  guides(colour = guide_colourbar(title.position="bottom", title.hjust = 0.5)) -> vwc_plot

png("figs/FigS6_vwc_rxn.png", height = 8, width = 8, res = 300, units = "in")
vwc_plot
dev.off()

## Fitness subplots (temp) ####

# Need to line up genotype indexing
genotype_ids_f <- as.numeric(as.factor(cg_model_fecun$genotype))
genotype_ids_s <- as.numeric(as.factor(cg_model$genotype))

tibble(genotype = factor(genotype_ids_f),
       genotype_id = cg_model_fecun$genotype) %>% 
  distinct() %>% 
  arrange(genotype) %>% 
  print(n = Inf)

tibble(genotype = factor(genotype_ids_s),
       genotype_id = cg_model$genotype) %>% 
  distinct() %>% 
  arrange(genotype) %>% 
  print(n = Inf)

# Ok these should be the same and will line up
ln_fitness_genotypes <- list()

for (i in 1:length(unique(genotype_ids_f))){
  ln_fitness_genotypes[[i]] <- log(exp(genotype_samples_f[[i]])*plogis(genotype_samples_s[[i]]))
}

result <- as.data.frame(sapply(ln_fitness_genotypes, function(df) rowMeans(df)))

names(result) <- paste("g", 1:96, sep = "_")
result %>% 
  mutate(temp = seq_temp) %>% 
  gather(key = genotype, value = ln_fitness, g_1:g_96) %>% 
  mutate(genotype = as.factor(parse_number(genotype))) -> fitness_preds

# Get genotype numeric codes that were used in fitting the model and line them
# back up with true codes
tibble(genotype = factor(genotype_ids),
       genotype_id = cg_model$genotype) %>% 
  distinct() %>% 
  merge(fitness_preds) %>% 
  dplyr::select(genotype = genotype_id, temp, ln_fitness) ->fitness_preds


# Get in PC1 values by genotype
cg_model %>% 
  dplyr::select(genotype, PC1) %>% 
  distinct() %>% 
  merge(fitness_preds) -> fitness_preds

fitness_preds %>% 
  ggplot(aes(x = temp, y = ln_fitness, group = genotype, color = PC1)) +
  geom_line(linewidth = 0.5) +
  theme_classic(base_size = 16) +
  labs(y = "ln(fitness)",
       x = "soil temperature (°C)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  ylim(2.4,4.5)-> temp_fitness_rxn

temp_fitness <- fitness_preds

temp_fitness %>% 
  filter(temp == min(temp_fitness$temp) | temp == max(temp_fitness$temp)) %>% 
  group_by(genotype) %>% 
  slice_min(temp) -> low_temp

temp_fitness %>% 
  filter(temp == min(temp_fitness$temp) | temp == max(temp_fitness$temp)) %>% 
  group_by(genotype) %>% 
  slice_max(temp) %>% 
  cbind(fitness2 = low_temp$ln_fitness) %>% 
  mutate(diff_fit = ln_fitness - fitness2) %>% 
  mutate(slope = diff_fit / (max(temp_fitness$temp) - min(temp_fitness$temp))) -> slope_fitness

coef(summary(lm(slope ~ PC1, data = slope_fitness)))[2,1] -> slope_fit
confint(lm(slope ~ PC1, data = slope_fitness))[2,] -> conf_int_fit

slope_fitness %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) +
  geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "ln(fitness)/", Delta, "temp."))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_hline(aes(yintercept = 0), color = "gray47", linewidth = 1.2) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed")+
  ylim(-0.07, 0.12) +
  # Adding annotation with slope and 95% CI
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf(
      "italic(beta)==%.3f*' [%.3f, %.3f]'",
      slope_fit, conf_int_fit[1], conf_int_fit[2]
    ),
    parse = TRUE,
    hjust = 1.1, vjust = 1.1,
    size = beta_size
  )-> temp_fitness_b


## Fitness subplot (vwc) ####
# Ok these should be the same and will line up
ln_fitness_genotypes <- list()

for (i in 1:length(unique(genotype_ids_f))){
  ln_fitness_genotypes[[i]] <- log(exp(genotype_samples_f[[i]])*plogis(genotype_samples_s[[i]]))
}

result <- as.data.frame(sapply(ln_fitness_genotypes, function(df) rowMeans(df)))

names(result) <- paste("g", 1:96, sep = "_")
result %>% 
  mutate(vwc = seq_vwc) %>% 
  gather(key = genotype, value = ln_fitness, g_1:g_96) %>% 
  mutate(genotype = as.factor(parse_number(genotype))) -> fitness_preds

# Get genotype numeric codes that were used in fitting the model and line them
# back up with true codes
tibble(genotype = factor(genotype_ids),
       genotype_id = cg_model$genotype) %>% 
  distinct() %>% 
  merge(fitness_preds) %>% 
  dplyr::select(genotype = genotype_id, vwc, ln_fitness) ->fitness_preds


# Get in PC1 values by genotype
cg_model %>% 
  dplyr::select(genotype, PC1) %>% 
  distinct() %>% 
  merge(fitness_preds) -> fitness_preds

# Rename object
vwc_fitness <- fitness_preds

vwc_fitness %>% 
  ggplot(aes(x = vwc, y = ln_fitness, group = genotype, color = PC1)) +
  geom_line(linewidth = 0.5) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "ln(fitness)",
       x = "soil moisture (vwc)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  ylim(2.4,4.5)-> vwc_fitness_rxn

vwc_fitness %>% 
  filter(vwc == min(vwc_fitness$vwc) | vwc == max(vwc_fitness$vwc)) %>% 
  group_by(genotype) %>% 
  slice_min(vwc) -> low_vwc_fecun

vwc_fitness %>% 
  filter(vwc == min(vwc_fitness$vwc) | vwc == max(vwc_fitness$vwc)) %>% 
  group_by(genotype) %>% 
  slice_max(vwc) %>% 
  cbind(fitness2 = low_vwc_fecun$ln_fitness) %>% 
  mutate(diff_fit = ln_fitness - fitness2) %>% 
  mutate(slope = diff_fit / (max(vwc_fitness$vwc) - min(vwc_fitness$vwc))) -> slope_fitness

coef(summary(lm(slope ~ PC1, data = slope_fitness)))[2,1] -> slope_fit
confint(lm(slope ~ PC1, data = slope_fitness))[2,] -> conf_int_fit

slope_fitness %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) +
  geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "ln(fitness)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_hline(aes(yintercept = 0), color = "gray47", linewidth = 1.2) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed")+
  # Adding annotation with slope and 95% CI
  annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf(
      "italic(beta)==%.3f*' [%.3f, %.3f]'",
      slope_fit, conf_int_fit[1], conf_int_fit[2]
    ),
    parse = TRUE,
    hjust = 1.1, vjust = 1.1,
    size = beta_size
  ) -> vwc_fitness_b



## Generate Figure 3 ####

temp_fitness_rxn + temp_fitness_b +
  vwc_fitness_rxn + vwc_fitness_b +
  plot_layout(guides = "collect", axis_titles = "collect_x") +
  plot_annotation(tag_levels = "a") &
  theme_classic(base_size = 14) &
  theme(legend.position = "bottom",
        legend.key.width=unit(1.5,"cm")) & 
  labs(color = "PC 1 (hot & dry → cool & wet)") &
  guides(colour = guide_colourbar(title.position="bottom", title.hjust = 0.5)) -> combined_fitness

png("figs/Fig3_fitness_combined.png", width = 8, height = 8, res = 300, units = "in")
combined_fitness
dev.off()
