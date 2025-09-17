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

# Get model draws
draws <- fit$draws(format = "df")
draws_s <- fit_s$draws(format = "df")

colors <- c("#d8b365", "#5ab4ac")

# Figure out which interactions are significant
quantile(draws$`beta[6]`, c(0.025, 0.975))# F: temp x nb = no
quantile(draws$`beta[7]`, c(0.025, 0.975))# F: vwc x nb = no
quantile(draws_s$`beta[6]`, c(0.025, 0.975))# S: temp x nb = yes
quantile(draws_s$`beta[7]`, c(0.025, 0.975))# S: vwc x nb = yes

## Subplot for vwc survival ####
cg_model$vwc_avg_sc <- scale(cg_model$vwc_avg)[,1]

seq_sc_vwc <- seq(min(cg_model$vwc_avg_sc),
                        max(cg_model$vwc_avg_sc),
                        length.out = 100)
seq_vwc<- seq(min(cg_model$vwc_avg),
                     max(cg_model$vwc_avg),
                     length.out = 100)

dens1_vs <- matrix(NA, nrow = 100, ncol = nrow(draws))
dens2_vs <- matrix(NA, nrow = 100, ncol = nrow(draws))

for (i in 1:100){
  dens1_vs[i,] <- draws_s$alpha +
    draws_s$`beta[1]` * -1 +
    draws_s$`beta[3]` * seq_sc_vwc[i] +
    draws_s$`beta[7]` * seq_sc_vwc[i] *-1 
  
  dens2_vs[i,] <- draws_s$alpha +
    draws_s$`beta[1]` +
    draws_s$`beta[3]` * seq_sc_vwc[i] +
    draws_s$`beta[7]` * seq_sc_vwc[i] 
}

plogis(apply(dens1_vs, 1, quantile, c(0.025, 0.5, 0.975))) -> high_dens_pred
plogis(apply(dens2_vs, 1, quantile, c(0.025, 0.5, 0.975))) -> low_dens_pred

# Calculate x value at peak of parabola
tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       survived = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("low density", "high density"), each = 100),
       vwc = c(seq_vwc, seq_vwc)) %>% 
  ggplot(aes(x = vwc, y = survived, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "P(survival)",
       x = "soil moisture (vwc)",
       fill = "",
       color = "") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylim(0.05, 0.95)-> survival_vwc

## Subplot for temp survival ####
cg_model$temp_surv_sc <- scale(cg_model$temp_surv)[,1]

seq_sc_temp <- seq(min(cg_model$temp_surv_sc),
                  max(cg_model$temp_surv_sc),
                  length.out = 100)
seq_temp<- seq(min(cg_model$temp_surv),
              max(cg_model$temp_surv),
              length.out = 100)

dens1_ts <- matrix(NA, nrow = 100, ncol = nrow(draws))
dens2_ts <- matrix(NA, nrow = 100, ncol = nrow(draws))

for (i in 1:100){
  dens1_ts[i,] <- draws_s$alpha +
    draws_s$`beta[1]` * -1 +
    draws_s$`beta[2]` * seq_sc_temp[i] +
    draws_s$`beta[6]` * seq_sc_temp[i] *-1 
  
  dens2_ts[i,] <- draws_s$alpha +
    draws_s$`beta[1]` +
    draws_s$`beta[2]` * seq_sc_temp[i] +
    draws_s$`beta[6]` * seq_sc_temp[i] 
}

plogis(apply(dens1_ts, 1, quantile, c(0.025, 0.5, 0.975))) -> high_dens_pred
plogis(apply(dens2_ts, 1, quantile, c(0.025, 0.5, 0.975))) -> low_dens_pred

tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       survived = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("low density", "high density"), each = 100),
       temp = c(seq_temp, seq_temp)) %>% 
  ggplot(aes(x = temp, y = survived, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "P(survival)",
       x = "soil temperature (°C)",
       fill = "",
       color = "") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)+
  ylim(0.05, 0.95)-> survival_temp

## Subplot for vwc fecundity ####

dens1_vf <- matrix(NA, nrow = 100, ncol = nrow(draws))
dens2_vf <- matrix(NA, nrow = 100, ncol = nrow(draws))

# Make data set of just plants that survived and reproduced
cg_model$survived <- ifelse(cg_model$seed_count > 0, 1, 0)

# Make subset of data that survived
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

# Figure out mean value of density (scaled) for each treatment
cg_model_fecun$new_neighbors_sc <- scale(cg_model_fecun$new_neighbors)
cg_model_fecun %>% 
  group_by(density) %>% 
  summarize(mean_dens = mean(new_neighbors_sc)) -> means_by_trt

cg_model_fecun$vwc_avg_sc <- scale(cg_model_fecun$vwc_avg)
seq_sc_vwc <- seq(min(cg_model_fecun$vwc_avg_sc),
                  max(cg_model_fecun$vwc_avg_sc),
                  length.out = 100)
seq_vwc<- seq(min(cg_model_fecun$vwc_avg),
              max(cg_model_fecun$vwc_avg),
              length.out = 100)

for (i in 1:100){
  dens1_vf[i,] <- draws$alpha +
    draws$`beta[1]` * means_by_trt %>% filter(density == "hi") %>% pull(mean_dens) +
    draws$`beta[3]` * seq_sc_vwc[i] +
    draws$`beta[7]` * seq_sc_vwc[i] * means_by_trt %>% filter(density == "hi") %>% pull(mean_dens)
  
  dens2_vf[i,] <- draws$alpha +
    draws$`beta[1]` * means_by_trt %>% filter(density == "lo") %>% pull(mean_dens) +
    draws$`beta[3]` * seq_sc_vwc[i] +
    draws$`beta[7]` * seq_sc_vwc[i] * means_by_trt %>% filter(density == "lo") %>% pull(mean_dens)
}

apply(dens1_vf, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_pred
apply(dens2_vf, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_pred

# Calculate x value at peak of parabola
tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       fecundity = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("low density", "high density"), each = 100),
       vwc = c(seq_vwc, seq_vwc)) %>% 
  ggplot(aes(x = vwc, y = fecundity, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "ln(fecundity)",
       x = "soil moisture (vwc)",
       fill = "",
       color = "") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) -> fecundity_vwc

## Subplot for temp fecundity ####
cg_model_fecun$temp_fecun_sc <- scale(cg_model_fecun$temp_fecun)[,1]

seq_sc_temp <- seq(min(cg_model_fecun$temp_fecun_sc),
                   max(cg_model_fecun$temp_fecun_sc),
                   length.out = 100)
seq_temp<- seq(min(cg_model_fecun$temp_fecun),
               max(cg_model_fecun$temp_fecun),
               length.out = 100)

dens1_tf <- matrix(NA, nrow = 100, ncol = nrow(draws))
dens2_tf <- matrix(NA, nrow = 100, ncol = nrow(draws))

for (i in 1:100){
  dens1_tf[i,] <- draws$alpha +
    draws$`beta[1]` * means_by_trt %>% filter(density == "hi") %>% pull(mean_dens) +
    draws$`beta[2]` * seq_sc_temp[i] +
    draws$`beta[6]` * seq_sc_temp[i] * means_by_trt %>% filter(density == "hi") %>% pull(mean_dens) 
  
  dens2_tf[i,] <- draws$alpha +
    draws$`beta[1]` * means_by_trt %>% filter(density == "lo") %>% pull(mean_dens) +
    draws$`beta[2]` * seq_sc_temp[i] +
    draws$`beta[6]` * seq_sc_temp[i] * means_by_trt %>% filter(density == "lo") %>% pull(mean_dens)
}

apply(dens1_tf, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_pred
apply(dens2_tf, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_pred

# Calculate x value at peak of parabola
tibble(lower = c(low_dens_pred[1,], high_dens_pred[1,]),
       upper = c(low_dens_pred[3,], high_dens_pred[3,]),
       fecundity = c(low_dens_pred[2,], high_dens_pred[2,]),
       density = rep(c("low density", "high density"), each = 100),
       temp = c(seq_temp, seq_temp)) %>% 
  ggplot(aes(x = temp, y = fecundity, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "ln(fecundity)",
       x = "soil temperature (°C)",
       fill = "",
       color = "") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) -> fecundity_temp

## Fitness for soil temperature ####
fitness_t_hi <- log(plogis(dens1_ts) * exp(dens1_tf))
fitness_t_lo <- log(plogis(dens2_ts) * exp(dens2_tf))

apply(fitness_t_hi, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_t_fit
apply(fitness_t_lo, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_t_fit

tibble(lower = c(low_dens_t_fit[1,], high_dens_t_fit[1,]),
       upper = c(low_dens_t_fit[3,], high_dens_t_fit[3,]),
       fitness = c(low_dens_t_fit[2,], high_dens_t_fit[2,]),
       density = rep(c("low density", "high density"), each = 100),
       temp = c(seq_temp, seq_temp)) %>% 
  ggplot(aes(x = temp, y = fitness, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "ln(fitness)",
       x = "soil temperature (°C)",
       fill = "",
       color = "") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) -> fitness_temp

## Fitness for soil moisture ###
fitness_v_hi <- log(plogis(dens1_vs) * exp(dens1_vf))
fitness_v_lo <- log(plogis(dens2_vs) * exp(dens2_vf))

apply(fitness_v_hi, 1, quantile, c(0.025, 0.5, 0.975)) -> high_dens_v_fit
apply(fitness_v_lo, 1, quantile, c(0.025, 0.5, 0.975)) -> low_dens_v_fit

tibble(lower = c(low_dens_v_fit[1,], high_dens_v_fit[1,]),
       upper = c(low_dens_v_fit[3,], high_dens_v_fit[3,]),
       fitness = c(low_dens_v_fit[2,], high_dens_v_fit[2,]),
       density = rep(c("low density", "high density"), each = 100),
       vwc = c(seq_vwc, seq_vwc)) %>% 
  ggplot(aes(x = vwc, y = fitness, color = density, group = density)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density), alpha = 0.2, color = NA) +
  theme_classic(base_size = 16) +
  labs(y = "ln(fitness)",
       x = "soil moisture (vwc)",
       fill = "",
       color = "") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) -> fitness_vwc


## Bring plots together ####

png("figs/Fig7_dens_climate.png", height = 8.5, width = 7.5, res = 300, units = "in")
(survival_temp & theme(plot.tag.position = c(0.02,0.85))) +
(survival_vwc & theme(plot.tag.position = c(0.02,0.85))) +
(fecundity_temp & theme(plot.tag.position = c(0.02,1))) +
(fecundity_vwc & theme(plot.tag.position = c(0.02,1))) +
  (fitness_temp & theme(plot.tag.position = c(0.02,1))) +
  (fitness_vwc & theme(plot.tag.position = c(0.02,1)))+
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", axis_titles = "collect", nrow = 3) &
  theme(legend.position = "top", legend.background=element_blank())
dev.off()
