## Preliminaries ####
# Load libraries
library(tidyverse); library(here); library(stringi)
library(prism); library(raster); library(geosphere); library(sf);
library(dismo); library(factoextra); library(ggfortify);
library(daymetr); library(lubridate); library(scales)

## Make predictions with model and prediction model matrices ####

# Read in model output
fit <- as_cmdstan_fit(files = c("outputs/demo_model_fecun-202506031125-1-391f4a.csv",
                                "outputs/demo_model_fecun-202506031125-2-391f4a.csv",
                                "outputs/demo_model_fecun-202506031125-3-391f4a.csv"))

fit_s <- as_cmdstan_fit(files = c("outputs/demo_model_surv_noncenter-202506031158-1-25d6a2.csv",
                                  "outputs/demo_model_surv_noncenter-202506031158-2-25d6a2.csv",
                                  "outputs/demo_model_surv_noncenter-202506031158-3-25d6a2.csv"))

# Collect intercept and temperature regression coefficient
draws <- fit$draws(format = "df")
# Global intercept
alpha <- mean(draws$alpha)
# Regression coefficient for temperature
beta_2 <- mean(draws$`beta[2]`)
# Regression coefficient for climate mismatch
beta_4 <- mean(draws$`beta[4]`)
# Regression coefficient for climate mismatch ^2
beta_5 <- mean(draws$`beta[5]`)
# Genotype intercepts
intercepts <- colMeans(draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix())
# Genotype x temp slopes
slopes <- colMeans(draws %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",2]")) %>% as.matrix())

# Read in matrix of covariates for now and the future
pred_matrix_source <- read_csv("outputs/predict_source_mm.csv")
pred_matrix_future <- read_csv("outputs/predict_future_mm.csv")

pred_now <- NULL
pred_future <- NULL

for (i in 1:nrow(pred_matrix_future)){
  pred_now[i] <- alpha + intercepts[i] + (beta_2 + slopes[i]) * pred_matrix_source$temp_fecun_sc[i]
  pred_future[i] <- alpha + intercepts[i] + (beta_2 + slopes[i]) * pred_matrix_future$temp_fecun_sc[i] +
    beta_4 * pred_matrix_future$clim_dist_sc_f[i] + beta_5 * pred_matrix_future$clim_dist_sc2_f[i]
}

## Repeat for survival ####
# Collect intercept and temperature regression coefficient
draws_s <- fit_s$draws(format = "df")
# Global intercept
alpha_s <- mean(draws_s$alpha)
# Regression coefficient for temperature
beta_2s <- mean(draws_s$`beta[2]`)
# Regression coefficient for climate mismatch
beta_4s <- mean(draws_s$`beta[4]`)
# Regression coefficient for climate mismatch ^2
beta_5s <- mean(draws_s$`beta[5]`)
# Genotype intercepts
intercepts_s <- colMeans(draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_intercepts")) %>% as.matrix())
# Genotype x temp slopes
slopes_s <- colMeans(draws_s %>% as_tibble() %>% dplyr::select(starts_with("genotype_slopes") & ends_with(",2]")) %>% as.matrix())

pred_now_s <- NULL
pred_future_s <- NULL

for (i in 1:nrow(pred_matrix_future)){
  pred_now_s[i] <- alpha_s + intercepts_s[i] + (beta_2s + slopes_s[i]) * pred_matrix_source$temp_surv_sc[i]
  pred_future_s[i] <- alpha_s + intercepts_s[i] + (beta_2s + slopes_s[i]) * pred_matrix_future$temp_surv_sc[i] +
    beta_4s * pred_matrix_future$clim_dist_sc_s[i] + beta_5s * pred_matrix_future$clim_dist_sc2_s[i]
}

predict_fitness_now <- exp(pred_now) * plogis(pred_now_s)
predict_fitness_future <- exp(pred_future) * plogis(pred_future_s)

# Determine color based on direction of change
arrow_colors <- ifelse(predict_fitness_future > predict_fitness_now, "red", "blue")

# Plotting data matrix
plot_data <- data.frame(
  soil_temp_now = pred_matrix_source$soil_temp_fecun,
  soil_temp_future = pred_matrix_source$soil_temp_fecun,
  fitness_now = log(predict_fitness_now),
  fitness_future = log(predict_fitness_future)
)

# Compute change and scaled magnitude
plot_data$fitness_change <- plot_data$fitness_future - plot_data$fitness_now
plot_data$fitness_change_scaled <- plot_data$fitness_change / max(abs(plot_data$fitness_change), na.rm = TRUE)

# Determine outline color: red for increase, blue for decrease
plot_data$outline_color <- ifelse(plot_data$fitness_change > 0, "red", "blue")

# Plot
ggplot() +
  # 1. Outline arrows (fixed color per row, NOT mapped via aes)
  geom_segment(data = plot_data,
               aes(x = soil_temp_now, y = fitness_now,
                   xend = soil_temp_future, yend = fitness_future),
               color = plot_data$outline_color,  # Set manually
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 2, lineend = "round") +
  
  # 2. Inner arrows with gradient fill (mapped to fitness_change_scaled)
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
