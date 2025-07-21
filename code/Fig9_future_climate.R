## Preliminaries ####
# Load libraries
library(tidyverse); library(here); library(stringi)
library(prism); library(raster); library(geosphere); library(sf);
library(dismo); library(factoextra); library(ggfortify);
library(daymetr); library(lubridate); library(scales)

## Make predictions with model and prediction model matrices ####

# Read in model output
fit <- as_cmdstan_fit(files = c("outputs/demo_model_fecun-202506091642-1-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-2-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-3-6d200d.csv"))

fit_s <- as_cmdstan_fit(files = c("outputs/demo_model_surv_noncenter-202506091738-1-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-2-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-3-51da52.csv"))

# Set seed
set.seed(4685)

# Collect intercept and temperature regression coefficient from fecundity model
draws <- fit$draws(format = "df")
# Global intercept
alpha <- draws$alpha
# Regression coefficient for temperature
beta_2 <- draws$`beta[2]`
# Regression coefficient for climate mismatch
beta_4 <- draws$`beta[4]`
# Regression coefficient for climate mismatch ^2
beta_5 <- draws$`beta[5]`
# Dispersion parameter for NB
phi <- draws$phi
# Genotype intercepts
intercepts <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()
# Genotype x temp slopes
slopes <- draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",2]")) %>% as.matrix()

# Read in matrix of covariates for now and the future
pred_matrix_source <- read_csv("outputs/predict_source_mm.csv")
pred_matrix_future <- read_csv("outputs/predict_future_mm.csv")

# Set up empty matrices to hold output
pred_now <- matrix(NA, nrow = nrow(intercepts), ncol = ncol(intercepts))
pred_future <- matrix(NA, nrow = nrow(intercepts), ncol = ncol(intercepts))
post_pred_now <- matrix(NA, nrow = nrow(intercepts), ncol = ncol(intercepts))
post_pred_future <- matrix(NA, nrow = nrow(intercepts), ncol = ncol(intercepts))

# Loop through iterations
for (i in 1:nrow(pred_now)){
  # Loop through genotypes
  for (j in 1:ncol(pred_now)){
    pred_now[i,j] <- alpha[i] + intercepts[i,j] + (beta_2[i] + slopes[i,j]) * pred_matrix_source$temp_fecun_sc[j]
    
    post_pred_now[i,j] <- rnbinom(1, size = phi[i], mu = exp(pred_now[i,j]))
    
    pred_future[i,j] <- alpha[i] + intercepts[i,j] + (beta_2[i] + slopes[i,j]) * pred_matrix_future$temp_fecun_sc[j] +
      beta_4[i] * pred_matrix_future$clim_dist_sc_f[j] + beta_5[i] * pred_matrix_future$clim_dist_sc2_f[j]
    
    post_pred_future[i,j] <- rnbinom(1, size = phi[i], mu = exp(pred_future[i,j]))
  }
}

## Repeat for survival ####
# Collect intercept and temperature regression coefficient
draws_s <- fit_s$draws(format = "df")
# Global intercept
alpha_s <- draws_s$alpha
# Regression coefficient for temperature
beta_2s <- draws_s$`beta[2]`
# Regression coefficient for climate mismatch
beta_4s <- draws_s$`beta[4]`
# Regression coefficient for climate mismatch ^2
beta_5s <- draws_s$`beta[5]`
# Genotype intercepts
intercepts_s <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix()
# Genotype x temp slopes
slopes_s <- draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",2]")) %>% as.matrix()

pred_now_s <- matrix(NA, nrow = nrow(intercepts_s), ncol = ncol(intercepts_s))
pred_future_s <- matrix(NA, nrow = nrow(intercepts_s), ncol = ncol(intercepts_s))
post_pred_now_s <- matrix(NA, nrow = nrow(intercepts_s), ncol = ncol(intercepts_s))
post_pred_future_s <- matrix(NA, nrow = nrow(intercepts_s), ncol = ncol(intercepts_s))

for (i in 1:nrow(pred_now_s)){
  for (j in 1:ncol(pred_now_s)){
    pred_now_s[i,j] <- alpha_s[i] + intercepts_s[i,j] + (beta_2s[i] + slopes_s[i,j]) * pred_matrix_source$temp_surv_sc[j]
    
    post_pred_now_s[i,j] <- rbinom(1, size = 1, prob = plogis(pred_now_s[i,j]))
    
    pred_future_s[i,j] <- alpha_s[i] + intercepts_s[i,j] + (beta_2s[i] + slopes_s[i,j]) * pred_matrix_future$temp_surv_sc[j] +
      beta_4s[i] * pred_matrix_future$clim_dist_sc_s[j] + beta_5s[i] * pred_matrix_future$clim_dist_sc2_s[j]
    
    post_pred_future_s[i,j] <- rbinom(1, size = 1, prob = plogis(pred_future_s[i,j]))
  }
}

# Define prediction without observation noise
predict_fitness_now <- log((miscTools::colMedians(exp(pred_now) * plogis(pred_now_s)))) 
predict_fitness_future <- log((miscTools::colMedians(exp(pred_future) * plogis(pred_future_s)))) 

# Determine color based on direction of change
arrow_colors <- ifelse(predict_fitness_future > predict_fitness_now, "red", "blue")

# Plotting data matrix
plot_data <- data.frame(
  soil_temp_now = pred_matrix_source$soil_temp_fecun,
  soil_temp_future = pred_matrix_source$soil_temp_fecun,
  fitness_now = predict_fitness_now,
  fitness_future = predict_fitness_future
)

# Compute change and scaled magnitude
plot_data$fitness_change <- plot_data$fitness_future - plot_data$fitness_now
plot_data$fitness_change_scaled <- plot_data$fitness_change / max(abs(plot_data$fitness_change), na.rm = TRUE)

# Determine outline color: red for increase, blue for decrease
plot_data$outline_color <- ifelse(plot_data$fitness_change > 0, "red", "blue")

# Plot
ggplot() +
  geom_segment(data = plot_data,
               aes(x = soil_temp_now, y = fitness_now,
                   xend = soil_temp_future, yend = fitness_future),
               color = plot_data$outline_color,  # Set manually
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 2, lineend = "round") +
  
  geom_segment(data = plot_data,
               aes(x = soil_temp_now, y = fitness_now,
                   xend = soil_temp_future, yend = fitness_future,
                   color = fitness_change_scaled),
               arrow = arrow(length = unit(0.15, "cm")),
               linewidth = 0.7, lineend = "round") +
  
  scale_color_gradientn(colors = c("blue", "white", "red"),
                        values = rescale(c(-1, 0, 1)),
                        limits = c(-1, 1),
                        guide = guide_colorbar(title = "Δ ln(fitness)")) +
  
  labs(x = "soil temperature (°C)\n[30-year climate normal]", y = "ln(fitness)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") -> future_climate

png("figs/Fig9_future_climate.png", height = 5.5, width = 6.6, units = "in", res = 300)
future_climate  
dev.off()
