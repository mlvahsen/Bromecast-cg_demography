library(scales)
## Fecundity ####
# Make vector of temperature scaled
seq_sc_vwc <- seq(min(cg_model_fecun$vwc_avg_sc),
                   max(cg_model_fecun$vwc_avg_sc),
                   length.out = 10)
# Make same length vector of temperature unscaled
seq_vwc <- seq(min(cg_model_fecun$vwc_avg),
                max(cg_model_fecun$vwc_avg),
                length.out = 10)

# Extract draws from model fit -- fecundity
draws <- fit$draws(format = "df")

# Collect intercept and temperature regression coefficient
alpha <- draws$alpha
beta_3 <- draws$`beta[3]`

# Precompute the intercepts and slopes for all genotypes (outside of loops)
intercepts <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()

# We assume that genotype_slopes is a matrix with dimensions: (samples x 3*genotypes)
slopes <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",3]")) %>% as.matrix()

# Initialize the output matrix
genotype_store<- matrix(NA, nrow = length(seq_temp), ncol = ncol(intercepts))

# Vectorized computation of genotype stores
for (j in 1:ncol(intercepts)) {
  # Extract the relevant intercept and slopes for genotype j
  intercept_j <- intercepts[, j]  # Vector of intercepts for genotype j
  
  # Correct indexing to get the slopes for the 2nd predictor for genotype j
  slope_j <- slopes[, j]  # Slopes for the second predictor (2nd column) for genotype j
  
  # Precompute the entire matrix of `store_temp_gs` for each genotype in a vectorized manner
  store_vwc_gs <- matrix(NA, nrow = length(seq_temp), ncol = nrow(intercepts))
  
  for (i in 1:length(seq_temp)) {
    store_vwc_gs[i, ] <- alpha + intercept_j + (beta_3 + slope_j) * seq_sc_vwc[i]
  }
  
  # Apply colMeans to reduce the result
  genotype_store[, j] <- rowMeans(store_vwc_gs)
}

as_tibble(genotype_store) -> genotype_store_df
names(genotype_store_df) <- paste("g", 1:90, sep = "_")
genotype_store_df %>% 
  mutate(vwc_avg = seq_vwc) %>% 
  gather(key = genotype, value = ln_seed_count, g_1:g_90) %>% 
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
  geom_line(linewidth = 1) +
  #geom_point(data = cg_model_fecun) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "ln(seed count)",
       x = "soil water content") +
  scale_color_distiller(palette = "RdYlBu", direction = 1)  -> a_vwc

vwc_g_pred %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = ln_seed_count[which(vwc_avg == min(vwc_avg))],
          max = ln_seed_count[which(vwc_avg == max(vwc_avg))],
          slope = (max - min) / (max(vwc_avg) - min(vwc_avg))) -> slope_fecun

summary(lm(slope ~ PC1, data = slope_fecun))
confint(lm(slope ~ PC1, data = slope_fecun))


slope_fecun %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1\n(hot & dry → cool & wet)",
       y =expression(paste(Delta, "ln(fecundity)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed")-> b_vwc

## Survival ####
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

# Vectorized computation of genotype stores
for (j in 1:ncol(intercepts)) {
  # Extract the relevant intercept and slopes for genotype j
  intercept_j <- intercepts[, j]  # Vector of intercepts for genotype j
  
  # Correct indexing to get the slopes for the 2nd predictor for genotype j
  slope_j <- slopes[, j]  # Slopes for the second predictor (2nd column) for genotype j
  
  # Precompute the entire matrix of `store_temp_gs` for each genotype in a vectorized manner
  store_vwc_gs <- matrix(NA, nrow = length(seq_sc_vwc), ncol = nrow(intercepts))
  
  for (i in 1:length(seq_sc_temp)) {
    store_vwc_gs[i, ] <- alpha + intercept_j + (beta_3 + slope_j) * seq_sc_vwc[i]
  }
  
  # Apply colMeans to reduce the result
  genotype_stores[, j] <- rowMeans(store_vwc_gs)
}

as_tibble(genotype_stores) -> genotype_store_dfs
names(genotype_store_dfs) <- paste("g", 1:90, sep = "_")
genotype_store_dfs %>% 
  mutate(vwc_avg = seq_vwc) %>% 
  gather(key = genotype, value = logit_survival, g_1:g_90) %>% 
  mutate(genotype = as.factor(parse_number(genotype)),
         survival = plogis(logit_survival)) -> vwc_g_preds

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
  geom_line(size = 1) +
  #geom_point(data = cg_model_fecun) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "P(survival)",
       x = "soil water content") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) -> a_vwc_s

vwc_g_preds %>% 
  group_by(genotype, PC1) %>% 
  reframe(min = logit_survival[which(vwc_avg == min(vwc_avg))],
          max = logit_survival[which(vwc_avg == max(vwc_avg))],
          slope = (max - min) / (max(vwc_avg) - min(vwc_avg))) -> slope_survival
  
summary(lm(slope ~ PC1, data = slope_survival))
confint(lm(slope ~ PC1, data = slope_survival))

slope_survival %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1\n(hot & dry → cool & wet)",
       y =expression(paste(Delta, "P(survival)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed")-> b_vwc_s

## Fitness ####
vwc_g_preds %>%
  distinct() %>% 
  arrange(genotype, vwc_avg) -> vwc_g_preds

vwc_g_pred %>%
  distinct() %>% 
  arrange(genotype, vwc_avg) -> vwc_g_pred

vwc_g_preds %>%
  cbind(vwc_fecun = vwc_g_pred$vwc_avg,
        seed_count = exp(vwc_g_pred$ln_seed_count)) %>% 
  mutate(fitness = survival * seed_count) -> vwc_fitness

vwc_fitness %>% 
  ggplot(aes(x = vwc_avg, y = log(fitness), group = genotype, color = PC1)) +
  geom_line(size = 1.2) +
  theme_classic(base_size = 16) +
  #theme(legend.position = "none") +
  labs(y = "ln(fitness)",
       x = "soil water content") +
  scale_color_distiller(palette = "RdYlBu", direction = 1) -> vwc_fitness_rxn

vwc_fitness %>% 
  filter(vwc_avg == min(vwc_fitness$vwc_avg) | vwc_avg == max(vwc_fitness$vwc_avg)) %>% 
  group_by(genotype) %>% 
  slice_min(vwc_avg) -> low_vwc_fecun

vwc_fitness %>% 
  filter(vwc_avg == min(vwc_fitness$vwc_avg) | vwc_avg == max(vwc_fitness$vwc_avg)) %>% 
  group_by(genotype) %>% 
  slice_max(vwc_avg) %>% 
  cbind(fitness2 = low_vwc_fecun$fitness) %>% 
  mutate(diff_fit = log(fitness) - log(fitness2)) %>% 
  mutate(slope = diff_fit / (max(vwc_fitness$vwc_avg) - min(vwc_fitness$vwc_avg))) -> slope_fitness

summary(lm(slope ~ PC1, data = slope_fitness))
confint(lm(slope ~ PC1, data = slope_fitness))

slope_fitness
  ggplot(aes(x = PC1, y = slope, color = PC1)) +
  geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1\n(hot & dry → cool & wet)",
       y =expression(paste(Delta, "ln(fitness)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed") -> vwc_fitness_b


## Bring plots together ####
design <- "AAAABB
           CCCCDD
           EEEEFF"
a_vwc_s + b_vwc_s +
  a_vwc + b_vwc +
  vwc_fitness_rxn + vwc_fitness_b +
  plot_layout(guides = "collect", design = design, axis_titles = "collect") +
  plot_annotation(tag_levels = "a") -> vwc_plot

png("figs/Fig4_vwc.png", width = 9, height = 9, res = 300, units = "in")
vwc_plot
dev.off()

survival <- a_temp_s + b_temp_s + plot_layout(guides = "collect", widths = c(2,1))
fecundity <- a_temp + b_temp + plot_layout(guides = "collect", widths = c(2,1))
fitness <- temp_fitness_rxn + temp_fitness_b + plot_layout(guides = "collect", widths = c(2,1))



survival / fecundity / fitness + plot_layout(guides = "collect")



