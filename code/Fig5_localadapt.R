# This code creates Figure 5: local adaptation evidence from climate mismatch

## Preliminaries ####

# Load libraries
library(scales); library(tidyverse); library(cmdstanr);
library(patchwork); library(grid)

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

# Get draws from the models
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

store_climdists <- matrix(NA, nrow = 100, ncol = nrow(draws_s))

for (i in 1:100){
  store_climdists[i,] <- draws_s$alpha + draws_s$`beta[4]` * seq_sc_climdists[i] +
    draws_s$`beta[5]` * seq_sc_climdists[i]^2 
}

# Calculate x value at peak of parabola
x_peak_s <- mean(-draws_s$`beta[4]`/(2*draws_s$`beta[5]`))
center_s <- attr(scale(cg_model$clim_dist), "scaled:center")
scale_s <- attr(scale(cg_model$clim_dist), "scaled:scale")

x_peak_s_sc <- x_peak_s*scale_s + center_s

plogis(apply(store_climdists, 1, quantile, c(0.025, 0.5, 0.975))) -> sum_stat_climdists

tibble(lower = sum_stat_climdists[1,],
       upper = sum_stat_climdists[3,],
       survived = sum_stat_climdists[2,],
       clim_dist = seq_climdists) %>% 
  ggplot(aes(x = clim_dist, y = survived)) +
  geom_jitter(data = cg_model, alpha = 0.2, height = 0.05, size = 0.1, shape = 16, width = 0.0005) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.8, fill = "gray47") +
  geom_line(linewidth = 1, color = "black") +
  theme_classic(base_size = 16) +
  labs(y = "P(survival)",
       x = "climate mismatch") +
  geom_vline(aes(xintercept = 0), color = "gray27", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = x_peak_s_sc), color = "brown1", linewidth = 1, linetype = "dashed") -> climdist_plot_s

## Subplot for fecundity ####
# Make subset of data that survived
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

cg_model_fecun$clim_dist_sc <- scale(cg_model_fecun$clim_dist)[,1]

seq_sc_climdist <- seq(min(cg_model_fecun$clim_dist_sc),
                       max(cg_model_fecun$clim_dist_sc),
                       length.out = 100)
seq_climdist <- seq(min(cg_model_fecun$clim_dist),
                    max(cg_model_fecun$clim_dist),
                    length.out = 100)

store_climdist <- matrix(NA, nrow = 100, ncol = nrow(draws))

for (i in 1:100){
  store_climdist[i,] <- draws$alpha + draws$`beta[4]` * seq_sc_climdist[i] +
    draws$`beta[5]` * seq_sc_climdist[i]^2 
}

apply(store_climdist, 1, quantile, c(0.025, 0.5, 0.975)) -> sum_stat_climdist

cg_model_fecun$ln_seed_count <- log(cg_model_fecun$seed_count)

# Calculate x value at peak of parabola
x_peak <- mean(-draws$`beta[4]`/(2*draws$`beta[5]`))

center <- attr(scale(cg_model_fecun$clim_dist), "scaled:center")
scale <- attr(scale(cg_model_fecun$clim_dist), "scaled:scale")

x_peak_sc <- x_peak*scale + center

tibble(lower = sum_stat_climdist[1,],
       upper = sum_stat_climdist[3,],
       ln_seed_count = sum_stat_climdist[2,],
       clim_dist = seq_climdist) %>% 
  ggplot(aes(x = clim_dist, y = ln_seed_count)) +
  geom_point(data = cg_model_fecun %>% filter(seed_count > 0), alpha = 0.2, size = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.8, fill = "gray47") +
  geom_line(linewidth = 1, color = "black") +
  theme_classic(base_size = 16) +
  labs(y = "ln(seed count)",
       x = "climate mismatch") +
  ylim(0,10)  +
  geom_vline(aes(xintercept = 0), color = "gray27", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = x_peak_sc), color = "dodgerblue", linewidth = 1, linetype = "dashed") -> climdist_plot

## Fitness subplot ####
fitness <- log(exp(sum_stat_climdist) * sum_stat_climdists)

text_high <- grid::textGrob("current environment\nhotter & drier\nthan source environment",
                            gp=gpar(fontsize=11, fontface="bold", lineheight = 0.8, col = "brown1"))
text_low <- grid::textGrob("current environment\ncooler & wetter\nthan source environment",
                           gp=gpar(fontsize=11, fontface="bold", lineheight = 0.8, col = "dodgerblue"))

# Find x where fitness is maximized
x_max_fit <- seq_climdist[which(fitness[2,] == max(fitness[2,]))]

tibble(lower = fitness[1,],
       upper = fitness[3,],
       ln_fitness = fitness[2,],
       clim_dist = seq_climdist) %>% 
  ggplot(aes(x = clim_dist, y = ln_fitness)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "gray47") +
  geom_line(linewidth = 1.5, color = "black") +
  theme_classic(base_size = 16) +
  labs(y = "ln(fitness)",
       x = "climate mismatch") +
  geom_vline(aes(xintercept = 0), color = "gray27", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = x_max_fit), color = "brown1", linewidth = 1, linetype = "dashed") +
  annotation_custom(text_high,xmin=6,xmax=6,ymin=-2.5,ymax=3) + 
  annotation_custom(text_low,xmin=-8,xmax=-8,ymin=-2.5,ymax=3)+
  coord_cartesian(ylim=c(1,5), clip="off") +
  coord_cartesian(clip = "off") + 
  ylim(1,5) +
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 20,  # Right margin
                             b = 30,  # Bottom margin
                             l = 10)) -> climdist_plot_fitness

## Bring plots together ####
design <- "AA
           AA
           BC"

climdist_plot_fitness + climdist_plot_s + climdist_plot +
  plot_layout(guides = "collect", design = design) +
  plot_annotation(tag_levels = "a") -> local_adaptation

png("figs/Fig5_localadapt.png", height = 8.25, width = 7.2, res = 300, units = "in")
local_adaptation
dev.off()

# Find biggest change in fitness
exp(min(fitness[2,])) # 12.5285
exp(min(fitness[1,])) # 3.647208
exp(min(fitness[3,])) # 40.82151

exp(max(fitness[2,])) # 42.0003
exp(max(fitness[1,])) # 15.5958
exp(max(fitness[3,])) # 109.9494
