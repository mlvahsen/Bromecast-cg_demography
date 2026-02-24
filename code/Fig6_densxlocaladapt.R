# This code creates Figure 6g-i which shows the interactions between climate
# mismatch and density dependence. The code generates this portion of Figure 6
# which is then used in 'Figure6_all.R' to create the full Figure 6 plot. This
# code also generates Figure S7 which shows the differences in curvature and
# optimal between high and low density treatments for fitness.

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

# Set plotting colors
colors <- c("#d8b365", "#5ab4ac")

# Get model draws
draws <- fit$draws(format = "df")
draws_s <- fit_s$draws(format = "df")

## Subplot for survival ####
cg_model$clim_dist_sc <- scale(cg_model$clim_dist)[,1]

seq_sc_climdists <- seq(min(cg_model$clim_dist_sc),
                        max(cg_model$clim_dist_sc),
                        length.out = 1000)
seq_climdists <- seq(min(cg_model$clim_dist),
                     max(cg_model$clim_dist),
                     length.out = 1000)

dens_low_s <- matrix(NA, nrow = 1000, ncol = nrow(draws))
dens_high_s <- matrix(NA, nrow = 1000, ncol = nrow(draws))

for (i in 1:1000){
  dens_high_s[i,] <- draws_s$alpha +
    draws_s$`beta[1]` * -1 +
    draws_s$`beta[4]` * seq_sc_climdists[i] +
    draws_s$`beta[5]` * seq_sc_climdists[i]^2 +
    draws_s$`beta[8]` * seq_sc_climdists[i] *-1 +
    draws_s$`beta[9]` * seq_sc_climdists[i]^2 *-1
  
  dens_low_s[i,] <- draws_s$alpha +
    draws_s$`beta[1]` +
    draws_s$`beta[4]` * seq_sc_climdists[i] +
    draws_s$`beta[5]` * seq_sc_climdists[i]^2 +
    draws_s$`beta[8]` * seq_sc_climdists[i] +
    draws_s$`beta[9]` * seq_sc_climdists[i]^2
}

plogis(apply(dens_high_s, 1, quantile, c(0.025, 0.5, 0.975))) -> high_dens_pred
plogis(apply(dens_low_s, 1, quantile, c(0.025, 0.5, 0.975))) -> low_dens_pred

# Calculate x value at peak of parabola
center_s <- attr(scale(cg_model$clim_dist), "scaled:center")
scale_s <- attr(scale(cg_model$clim_dist), "scaled:scale")

x_peak_s_low <- mean(-(draws_s$`beta[4]`+draws_s$`beta[8]`)/(2*(draws_s$`beta[5]`+draws_s$`beta[9]`)))
x_peak_s_low_sc <- x_peak_s_low*scale_s + center_s
x_peak_s_high <- mean(-(draws_s$`beta[4]`-draws_s$`beta[8]`)/(2*(draws_s$`beta[5]`-draws_s$`beta[9]`)))
x_peak_s_high_sc <- x_peak_s_high*scale_s + center_s

tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       survived = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("low density", "high density"), each = 1000),
       clim_dist = c(seq_climdists, seq_climdists)) %>% 
  ggplot(aes(x = clim_dist, y = survived, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "P(survival)",
       x = "climate mismatch",
       fill = "",
       color = "") +
  geom_vline(aes(xintercept = 0), color = "gray27", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = x_peak_s_high_sc), color = colors[1], linewidth = 1) +
  geom_vline(aes(xintercept = x_peak_s_low_sc), color = colors[2], linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) -> survival_sub

## Subplot for fecundity ####
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

seq_sc_climdist <- seq(min(cg_model_fecun$clim_dist_sc),
                       max(cg_model_fecun$clim_dist_sc),
                       length.out = 1000)
seq_climdist <- seq(min(cg_model_fecun$clim_dist),
                    max(cg_model_fecun$clim_dist),
                    length.out = 1000)

dens_high <- matrix(NA, nrow = 1000, ncol = nrow(draws))
dens_low <- matrix(NA, nrow = 1000, ncol = nrow(draws))

# Two different parabolas will be the average densities for each original
# density treatment
cg_model_fecun$new_neighbors_sc <- scale(cg_model_fecun$new_neighbors)[,1]
low_dens <- cg_model_fecun %>% filter(density == "lo") %>% pull(new_neighbors_sc) %>% mean()
high_dens <- cg_model_fecun %>% filter(density == "hi") %>% pull(new_neighbors_sc) %>% mean()

for (i in 1:1000){
  dens_low[i,] <- draws$alpha +
    draws$`beta[1]` * low_dens +
    draws$`beta[4]` * seq_sc_climdist[i] +
    draws$`beta[5]` * seq_sc_climdist[i]^2 +
    draws$`beta[8]` * seq_sc_climdist[i] * low_dens +
    draws$`beta[9]` * seq_sc_climdist[i]^2 * low_dens
  
  dens_high[i,] <- draws$alpha +
    draws$`beta[1]` * high_dens +
    draws$`beta[4]` * seq_sc_climdist[i] +
    draws$`beta[5]` * seq_sc_climdist[i]^2 +
    draws$`beta[8]` * seq_sc_climdist[i] * high_dens +
    draws$`beta[9]` * seq_sc_climdist[i]^2 * high_dens
}

apply(dens_low, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_pred
apply(dens_high, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_pred

# Calculate x value at peak of parabola
center <- attr(scale(cg_model$clim_dist), "scaled:center")
scale <- attr(scale(cg_model$clim_dist), "scaled:scale")

x_peak_low <- mean(-(draws$`beta[4]`+draws$`beta[8]`)/(2*(draws$`beta[5]`+draws$`beta[9]`)))
x_peak_low_sc <- x_peak_low*scale + center
x_peak_high <- mean(-(draws$`beta[4]`-draws$`beta[8]`)/(2*(draws$`beta[5]`-draws$`beta[9]`)))
x_peak_high_sc <- x_peak_high*scale + center

text_high <- grid::textGrob("current environment\nhotter & drier\nthan source environment",
                            gp=grid::gpar(fontsize=9, fontface="bold", lineheight = 0.8, col = "brown1"))
text_low <- grid::textGrob("current environment\ncooler & wetter\nthan source environment",
                           gp=grid::gpar(fontsize=9, fontface="bold", lineheight = 0.8, col = "dodgerblue"))

tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       fecundity = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("low density", "high density"), each = 1000),
       clim_dist = c(seq_climdists, seq_climdists)) %>% 
  ggplot(aes(x = clim_dist, y = fecundity, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "ln(fecundity)",
       x = "climate mismatch",
       fill = "",
       color = "") +
  geom_vline(aes(xintercept = 0), color = "gray27", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = x_peak_high_sc), color = colors[1], linewidth = 1) +
  geom_vline(aes(xintercept = x_peak_low_sc), color = colors[2], linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)  -> fecundity_sub

## Fitness plot ####
fitness_pred_high <- log(plogis(dens_high_s) * exp(dens_high))
fitness_pred_low <- log(plogis(dens_low_s) * exp(dens_low))

apply(fitness_pred_high, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_pred_fit
apply(fitness_pred_low, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_pred_fit

tibble(lower = c(low_dens_pred_fit[1,], high_dens_pred_fit[1,]),
       upper = c(low_dens_pred_fit[3,], high_dens_pred_fit[3,]),
       fitness = c(low_dens_pred_fit[2,], high_dens_pred_fit[2,]),
       density = rep(c("low density", "high density"), each = 1000),
       clim_dist = c(seq_climdists, seq_climdists)) -> fit_data  

# Find the highest value of fitness
fit_data %>% 
  group_by(density) %>% 
  summarize(max = max(fitness)) -> pred_max
  
x_low_peak <- as.numeric(fit_data[which(fit_data$fitness == pred_max$max[2]), "clim_dist"])
x_high_peak <- as.numeric(fit_data[which(fit_data$fitness == pred_max$max[1]), "clim_dist"])

fit_data %>% 
  ggplot(aes(x = clim_dist, y = fitness, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "ln(fitness)",
       x = "climate mismatch",
       fill = "",
       color = "") +
  geom_vline(aes(xintercept = 0), color = "black", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = x_high_peak), color = colors[1], linewidth = 1) +
  geom_vline(aes(xintercept = x_low_peak), color = colors[2], linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme(plot.margin = margin(t = 20,  # Top margin
                             r = 30,  # Right margin
                             b = 20,  # Bottom margin
                             l = 10)) -> fitness_sub


## Generate climate mismatch x density subplots for Figure 6 ####
saveRDS(survival_sub, "outputs/survival_sub.rds")
saveRDS(fecundity_sub, "outputs/fecundity_sub.rds")
saveRDS(fitness_sub, "outputs/fitness_sub.rds")

## Generate Figure S7 ####
quad_slope_high <- NULL
quad_slope_low <- NULL
lin_slope_high <- NULL
lin_slope_low <- NULL

for (i in 1:ncol(fitness_pred_high)){
  quad_slope_high[i] <- as.numeric(coef(lm(fitness_pred_high[,i] ~ seq_climdists + I(seq_climdists^2)))[3])
  quad_slope_low[i] <- as.numeric(coef(lm(fitness_pred_low[,i] ~ seq_climdists + I(seq_climdists^2)))[3])
  lin_slope_high[i] <- as.numeric(coef(lm(fitness_pred_high[,i] ~ seq_climdists + I(seq_climdists^2)))[2])
  lin_slope_low[i] <- as.numeric(coef(lm(fitness_pred_low[,i] ~ seq_climdists + I(seq_climdists^2)))[2])
}

tibble(diff = quad_slope_high - quad_slope_low) %>% 
  ggplot(aes(x = diff)) +
  geom_density(fill = "gray") +
  xlab("difference in curvature\n(high density - low density)") +
  theme_classic(base_size = 16) +
  geom_vline(aes(xintercept = quantile(diff, 0.025)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(diff, 0.975)), linetype = "dashed") -> a_fit
  
quantile(quad_slope_high - quad_slope_low, c(0.025, 0.975))
# 2.5%         97.5% 
# -0.0007397001  0.0070221524 

opt_high <- -lin_slope_high / (2*quad_slope_high)
opt_low <- -lin_slope_low / (2*quad_slope_low)

tibble(opt = opt_high - opt_low) %>% 
  ggplot(aes(x = opt)) +
  geom_density(fill = "gray") +
  xlab("difference in x-value at maximum\n(high density - low density)") +
  theme_classic(base_size = 15) +
  geom_vline(aes(xintercept = quantile(opt, 0.025)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(opt, 0.975)), linetype = "dashed") -> b_fit

quantile(opt_high - opt_low, c(0.025, 0.975))
# 2.5%      97.5% 
# 0.03978238 2.57544456 

png("figs/FigS7_dens_adapt.png", height = 4, width = 7.5, units = "in", res = 300)
a_fit + b_fit
dev.off()
