library(scales)
## Fecundity ####
# Make vector of temperature scaled
seq_sc_temp <- seq(min(cg_model_fecun$temp_fecun_sc),
                   max(cg_model_fecun$temp_fecun_sc),
                   length.out = 10)
# Make same length vector of temperature unscaled
seq_temp <- seq(min(cg_model_fecun$temp_fecun),
                max(cg_model_fecun$temp_fecun),
                length.out = 10)

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
}

as_tibble(genotype_store) -> genotype_store_df
names(genotype_store_df) <- paste("g", 1:90, sep = "_")
genotype_store_df %>% 
  mutate(temp_fecun = seq_temp) %>% 
  gather(key = genotype, value = ln_seed_count, g_1:g_90) %>% 
  mutate(genotype = as.factor(parse_number(genotype))) -> temp_g_pred

# Get genotype numeric codes that were used in fitting the model and line them
# back up with true codes
genotype_id <- as.numeric(as.factor(cg_model_fecun$genotype))

tibble(genotype = factor(genotype_id),
       genotype_id = cg_model_fecun$genotype) %>% 
  distinct() %>% 
  merge(temp_g_pred) %>% 
  dplyr::select(genotype = genotype_id, temp_fecun, ln_seed_count) ->temp_g_pred

# Collect PC1 values for each genotype
cg_model %>% 
  dplyr::select(genotype, PC1) %>% 
  distinct() %>% 
  merge(temp_g_pred) -> temp_g_pred

temp_g_pred %>% 
  ggplot(aes(x = temp_fecun, y = ln_seed_count, group = genotype, color = PC1)) +
  geom_line(linewidth = 1) +
  #geom_point(data = cg_model_fecun) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "ln(seed count)",
       x = "temperature (°C)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1)  -> a_temp

temp_g_pred %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = ln_seed_count[which(temp_fecun == min(temp_fecun))],
          max = ln_seed_count[which(temp_fecun == max(temp_fecun))],
          slope = (max - min) / (max(temp_fecun) - min(temp_fecun))) -> slope_fecun

summary(lm(slope ~ PC1, data = slope_fecun))
confint(lm(slope ~ PC1, data = slope_fecun))

slope_fecun %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1\n(hot & dry → cool & wet)",
       y =expression(paste(Delta, "ln(fecundity)/", Delta, "temp."))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed")-> b_temp

## Survival ####
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
seq_sc_temps <- seq_sc_temp

# Precompute the intercepts and slopes for all genotypes (outside of loops)
intercepts <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()

# We assume that genotype_slopes is a matrix with dimensions: (samples x 3*genotypes)
slopes <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",2]")) %>% as.matrix()

# Initialize the output matrix
genotype_stores <- matrix(NA, nrow = length(seq_sc_temps), ncol = ncol(intercepts))

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
}

as_tibble(genotype_stores) -> genotype_store_dfs
names(genotype_store_dfs) <- paste("g", 1:90, sep = "_")
genotype_store_dfs %>% 
  mutate(temp_surv = seq_temps) %>% 
  gather(key = genotype, value = logit_survival, g_1:g_90) %>% 
  mutate(genotype = as.factor(parse_number(genotype)),
         survival = plogis(logit_survival)) -> temp_g_preds

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
  geom_line(size = 1) +
  #geom_point(data = cg_model_fecun) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "P(survival)",
       x = "temperature (°C)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) -> a_temp_s

temp_g_preds %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = logit_survival[which(temp_surv == min(temp_surv))],
          max = logit_survival[which(temp_surv == max(temp_surv))],
          slope = (max - min) / (max(temp_surv) - min(temp_surv))) ->slope_surv
  
summary(lm(slope ~ PC1, data = slope_surv))
confint(lm(slope ~ PC1, data = slope_surv))


slope_surv %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1\n(hot & dry → cool & wet)",
       y =expression(paste(Delta, "P(survival)/", Delta, "temp."))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed")-> b_temp_s

## Fitness ####
temp_g_preds %>%
  distinct() %>% 
  arrange(genotype, temp_surv) -> temp_g_preds

temp_g_pred %>%
  distinct() %>% 
  arrange(genotype, temp_fecun) -> temp_g_pred

temp_g_preds %>%
  cbind(temp_fecun = temp_g_pred$temp_fecun,
        seed_count = exp(temp_g_pred$ln_seed_count)) %>% 
  mutate(fitness = survival * seed_count) -> temp_fitness
  
temp_fitness %>% 
  ggplot(aes(x = temp_fecun, y = log(fitness), group = genotype, color = PC1)) +
  geom_line(size = 1.2) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "ln(fitness)",
       x = "temperature (°C)") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) -> temp_fitness_rxn

temp_fitness %>% 
  filter(temp_fecun == min(temp_fitness$temp_fecun) | temp_fecun == max(temp_fitness$temp_fecun)) %>% 
  group_by(genotype) %>% 
  slice_min(temp_fecun) -> low_temp_fecun

temp_fitness %>% 
  filter(temp_fecun == min(temp_fitness$temp_fecun) | temp_fecun == max(temp_fitness$temp_fecun)) %>% 
  group_by(genotype) %>% 
  slice_max(temp_fecun) %>% 
  cbind(fitness2 = low_temp_fecun$fitness) %>% 
  mutate(diff_fit = log(fitness) - log(fitness2)) %>% 
  mutate(slope = diff_fit / (max(temp_fitness$temp_fecun) - min(temp_fitness$temp_fecun))) -> slope_fitness

summary(lm(slope ~ PC1, data = slope_fitness))
confint(lm(slope ~ PC1, data = slope_fitness))

slope_fitness %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) +
  geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1\n(hot & dry → cool & wet)",
       y =expression(paste(Delta, "ln(fitness)/", Delta, "temp."))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed") -> temp_fitness_b



## Bring plots together ####
design <- "AAAABB
           CCCCDD
           EEEEFF"
a_temp_s + b_temp_s +
  a_temp + b_temp +
  temp_fitness_rxn + temp_fitness_b +
  plot_layout(guides = "collect", design = design, axis_titles = "collect") +
  plot_annotation(tag_levels = "a") -> temp_plot

png("figs/Fig3_temp.png", width = 9, height = 9, res = 300, units = "in")
temp_plot
dev.off()
