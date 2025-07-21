# This code creates Figure 3: G × moisture reaction norms and sensitivities

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
cg_model_fecun$temp_fecun_sc <- scale(cg_model_fecun$temp_fecun)[,1]
cg_model_fecun$vwc_avg_sc <- scale(sqrt(cg_model_fecun$vwc_avg))[,1]
cg_model_fecun$temp_surv_sc <- scale(sqrt(cg_model_fecun$temp_surv))[,1]
cg_model_fecun$genotype <- as.factor(cg_model_fecun$genotype)
cg_model_fecun$survived <- ifelse(cg_model_fecun$survived == "Y", 1, 0)
cg_model_fecun$site_year_gravel <- as.factor(paste(cg_model_fecun$site,
                                                   cg_model_fecun$year,
                                                   cg_model_fecun$albedo, sep = "_"))

# Make vector of soil moisture scaled
seq_sc_vwc <- seq(min(cg_model_fecun$vwc_avg_sc),
                   max(cg_model_fecun$vwc_avg_sc),
                   length.out = 10)
# Make same length vector of soil moisture unscaled
seq_vwc <- seq(min(cg_model_fecun$vwc_avg),
                max(cg_model_fecun$vwc_avg),
                length.out = 10)

## Fecundity subplot ####
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

summary(lm(slope ~ PC1, data = slope_fecun))
confint(lm(slope ~ PC1, data = slope_fecun))

# Set size for text annotation
beta_size <- 4

slope_fecun %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "ln(fecundity)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +
  # Adding annotation with slope and 95% CI
  annotate("text", x = Inf, y = Inf, label = {
    # Fit the linear model
    model <- lm(slope ~ PC1, data = slope_fecun)
    # Extract slope and 95% CI
    slope <- coef(model)[2]
    conf_int <- confint(model)[2, ]
    # Format the label with β symbol and CI
    bquote(italic(beta) == .(round(slope, 3)) ~ 
             .(paste("[", round(conf_int[1], 3), ",", format(round(conf_int[2], 3), nsmall = 3), "]", sep = "")))
  }, hjust = 1.1, vjust = 1.1, size = beta_size, color = "black") +
  scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(5, 14))-> b_vwc

## Survival subplot ####

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
  
summary(lm(slope ~ PC1, data = slope_survival))
confint(lm(slope ~ PC1, data = slope_survival))

slope_survival %>% 
  ggplot(aes(x = PC1, y = slope, color = PC1)) + geom_point(size = 2) +
  theme_classic(base_size = 16) +
  labs(x = "PC 1",
       y =expression(paste(Delta, "P(survival)/", Delta, "vwc"))) +
  scale_color_distiller(palette = "RdYlBu", direction = 1) +
  geom_smooth(method = "lm", se = T, color = "black", linetype = "dashed") +
  # Adding annotation with slope and 95% CI
  annotate("text", x = Inf, y = Inf, label = {
    # Fit the linear model
    model <- lm(slope ~ PC1, data = slope_survival)
    # Extract slope and 95% CI
    slope <- coef(model)[2]
    conf_int <- confint(model)[2, ]
    # Format the label with β symbol and CI
    bquote(italic(beta) == .(format(round(slope, 3), nsmall = 3)) ~ 
             .(paste("[", round(conf_int[1], 3), ",", format(round(conf_int[2], 3), nsmall = 3), "]", sep = "")))
  }, hjust = 1.1, vjust = 1.1, size = beta_size, color = "black")+
  scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(-20,-4))-> b_vwc_s

## Fitness subplot ####
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
  scale_color_distiller(palette = "RdYlBu", direction = 1) -> vwc_fitness_rxn

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

summary(lm(slope ~ PC1, data = slope_fitness))
confint(lm(slope ~ PC1, data = slope_fitness))

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
  annotate("text", x = Inf, y = Inf, label = {
    # Fit the linear model
    model <- lm(slope ~ PC1, data = slope_fitness)
    # Extract slope and 95% CI
    slope <- coef(model)[2]
    conf_int <- confint(model)[2, ]
    # Format the label with β symbol and CI
    bquote(italic(beta) == .(format(round(slope, 3), nsmall = 3)) ~ 
             .(paste("[", round(conf_int[1], 3), ",", format(round(conf_int[2], 3), nsmall = 3), "]", sep = "")))
  }, hjust = 1.1, vjust = 1.1, size = beta_size, color = "black") + 
  scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(-3,8)) -> vwc_fitness_b


## Bring plots together ####
design <- "AAAABB
           CCCCDD
           EEEEFF"
a_vwc_s + b_vwc_s +
  a_vwc + b_vwc +
  vwc_fitness_rxn + vwc_fitness_b +
  plot_layout(guides = "collect", design = design, axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        legend.key.width=unit(1.5,"cm")) & 
  labs(color = "PC 1 (hot & dry → cool & wet)") &
  guides(colour = guide_colourbar(title.position="bottom", title.hjust = 0.5))-> vwc_plot

png("figs/Fig4_vwc.png", width = 9, height = 10, res = 300, units = "in")
vwc_plot
dev.off()