# This code creates Figure 3: G × temperature reaction norms and sensitivities

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
                        length.out = 100)
seq_climdists <- seq(min(cg_model$clim_dist),
                     max(cg_model$clim_dist),
                     length.out = 100)

dens1_s <- matrix(NA, nrow = 100, ncol = nrow(draws))
dens2_s <- matrix(NA, nrow = 100, ncol = nrow(draws))

for (i in 1:100){
  dens1_s[i,] <- draws_s$alpha +
    draws_s$`beta[1]` * -1 +
    draws_s$`beta[4]` * seq_sc_climdists[i] +
    draws_s$`beta[5]` * seq_sc_climdists[i]^2 +
    draws_s$`beta[8]` * seq_sc_climdists[i] *-1 +
    draws_s$`beta[9]` * seq_sc_climdists[i]^2 *-1
  
  dens2_s[i,] <- draws_s$alpha +
    draws_s$`beta[1]` +
    draws_s$`beta[4]` * seq_sc_climdists[i] +
    draws_s$`beta[5]` * seq_sc_climdists[i]^2 +
    draws_s$`beta[8]` * seq_sc_climdists[i] +
    draws_s$`beta[9]` * seq_sc_climdists[i]^2
}

plogis(apply(dens1_s, 1, quantile, c(0.025, 0.5, 0.975))) -> low_dens_pred
plogis(apply(dens2_s, 1, quantile, c(0.025, 0.5, 0.975))) -> high_dens_pred

# Calculate x value at peak of parabola
x_peak_s_low <- mean(-(draws_s$`beta[4]`+draws_s$`beta[8]`)/(2*(draws_s$`beta[5]`+draws_s$`beta[9]`)))
center_s <- attr(scale(cg_model$clim_dist), "scaled:center")
scale_s <- attr(scale(cg_model$clim_dist), "scaled:scale")

x_peak_s_low_sc <- x_peak_s_low*scale_s + center_s
x_peak_s_high <- median(-(draws_s$`beta[4]`-draws_s$`beta[8]`)/(2*(draws_s$`beta[5]`-draws_s$`beta[9]`)))
x_peak_s_high_sc <- x_peak_s_high*scale_s + center_s

tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       survived = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("high density", "low density"), each = 100),
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
                       length.out = 100)
seq_climdist <- seq(min(cg_model_fecun$clim_dist),
                    max(cg_model_fecun$clim_dist),
                    length.out = 100)

dens1 <- matrix(NA, nrow = 100, ncol = nrow(draws))
dens2 <- matrix(NA, nrow = 100, ncol = nrow(draws))

# Two different parabolas will be the average densities for each original
# density treatment
cg_model_fecun$new_neighbors_sc <- scale(cg_model_fecun$new_neighbors)[,1]
min_dens <- cg_model_fecun %>% filter(density == "lo") %>% pull(new_neighbors_sc) %>% mean()
max_dens <- cg_model_fecun %>% filter(density == "hi") %>% pull(new_neighbors_sc) %>% mean()

for (i in 1:100){
  dens1[i,] <- draws$alpha +
    draws$`beta[1]` * min_dens +
    draws$`beta[4]` * seq_sc_climdists[i] +
    draws$`beta[5]` * seq_sc_climdists[i]^2 +
    draws$`beta[8]` * seq_sc_climdists[i] * min_dens +
    draws$`beta[9]` * seq_sc_climdists[i]^2 * min_dens
  
  dens2[i,] <- draws$alpha +
    draws$`beta[1]` * max_dens +
    draws$`beta[4]` * seq_sc_climdists[i] +
    draws$`beta[5]` * seq_sc_climdists[i]^2 +
    draws$`beta[8]` * seq_sc_climdists[i] * max_dens +
    draws$`beta[9]` * seq_sc_climdists[i]^2 * max_dens
}

apply(dens1, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_pred
apply(dens2, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_pred

# Calculate x value at peak of parabola
x_peak_low <- mean(-(draws$`beta[4]`+draws$`beta[8]`)/(2*(draws$`beta[5]`+draws$`beta[9]`)))
center <- attr(scale(cg_model$clim_dist), "scaled:center")
scale <- attr(scale(cg_model$clim_dist), "scaled:scale")

x_peak_low_sc <- x_peak_low*scale + center
x_peak_high <- median(-(draws$`beta[4]`-draws$`beta[8]`)/(2*(draws$`beta[5]`-draws$`beta[9]`)))
x_peak_high_sc <- x_peak_high*scale + center

text_high <- grid::textGrob("current environment\nhotter & drier\nthan source environment",
                            gp=grid::gpar(fontsize=9, fontface="bold", lineheight = 0.8, col = "brown1"))
text_low <- grid::textGrob("current environment\ncooler & wetter\nthan source environment",
                           gp=grid::gpar(fontsize=9, fontface="bold", lineheight = 0.8, col = "dodgerblue"))

tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       fecundity = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("low density", "high density"), each = 100),
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
fitness_pred_1 <- log(plogis(dens1_s) * exp(dens2))
fitness_pred_2 <- log(plogis(dens2_s) * exp(dens1))

apply(fitness_pred_1, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_pred_fit
apply(fitness_pred_2, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_pred_fit

tibble(lower = c(low_dens_pred_fit[1,], high_dens_pred_fit[1,]),
       upper = c(low_dens_pred_fit[3,], high_dens_pred_fit[3,]),
       fitness = c(low_dens_pred_fit[2,], high_dens_pred_fit[2,]),
       density = rep(c("low density", "high density"), each = 100),
       clim_dist = c(seq_climdists, seq_climdists)) -> fit_data  

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
  annotation_custom(text_high,xmin=7,xmax=7,ymin=0,ymax=0.3) + 
  annotation_custom(text_low,xmin=-9,xmax=-9,ymin=0,ymax=0.3)+
  coord_cartesian(ylim=c(1,5), clip="off") +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(t = 20,  # Top margin
                             r = 30,  # Right margin
                             b = 20,  # Bottom margin
                             l = 10)) -> fitness_sub


## Full plot ####
png("figs/Fig7_densxlocaladapt.png", height = 9.5, width = 5.5, res = 300, units = "in")
(survival_sub & theme(plot.tag.position = c(0, 0.9))) / (fecundity_sub & theme(plot.tag.position = c(0.0, 1.05))) /
  (fitness_sub & theme(plot.tag.position = c(0.0, 1.05)))+
  plot_annotation(tag_levels = "a") +
  plot_layout(axis_titles = "collect", guides = "collect") &
  theme(legend.position = "top")
dev.off()

quad_slope_1 <- NULL
quad_slope_2 <- NULL
lin_slope_1 <- NULL
lin_slope_2 <- NULL

for (i in 1:ncol(fitness_pred_1)){
  quad_slope_1[i] <- as.numeric(coef(lm(log(fitness_pred_1[,i]) ~ seq_climdists + I(seq_climdists^2)))[3])
  quad_slope_2[i] <- as.numeric(coef(lm(log(fitness_pred_2[,i]) ~ seq_climdists + I(seq_climdists^2)))[3])
  lin_slope_1[i] <- as.numeric(coef(lm(log(fitness_pred_1[,i]) ~ seq_climdists + I(seq_climdists^2)))[2])
  lin_slope_2[i] <- as.numeric(coef(lm(log(fitness_pred_2[,i]) ~ seq_climdists + I(seq_climdists^2)))[2])
}

tibble(diff = quad_slope_1 - quad_slope_2) %>% 
  ggplot(aes(x = diff)) +
  geom_density(fill = "gray") +
  xlab("difference in curvature\n(high density - low density)") +
  theme_classic(base_size = 16) +
  geom_vline(aes(xintercept = quantile(diff, 0.025)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(diff, 0.975)), linetype = "dashed") -> a_fit
  
quantile(quad_slope_1 - quad_slope_2, c(0.025, 0.975))

opt1 <- -lin_slope_1 / (2*quad_slope_1)
opt2 <- -lin_slope_2 / (2*quad_slope_2)

tibble(opt = opt1 - opt2) %>% 
  ggplot(aes(x = opt)) +
  geom_density(fill = "gray") +
  xlab("difference in optimum\n(high density - low density)") +
  theme_classic(base_size = 16) +
  geom_vline(aes(xintercept = quantile(opt, 0.025)), linetype = "dashed") +
  geom_vline(aes(xintercept = quantile(opt, 0.975)), linetype = "dashed") -> b_fit

quantile(opt1-opt2, c(0.025, 0.975))

png("figs/FigS4_dens_adapt.png", height = 4, width = 7.5, units = "in", res = 300)
a_fit + b_fit
dev.off()