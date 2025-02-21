# Figure 5 local adaptation

## Subplot for survival ####
cg_model$clim_dist_sc <- scale(cg_model$clim_dist)[,1]

seq_sc_climdists <- seq(min(cg_model$clim_dist_sc),
                        max(cg_model$clim_dist_sc),
                        length.out = 10)
seq_climdists <- seq(min(cg_model$clim_dist),
                     max(cg_model$clim_dist),
                     length.out = 10)

store_climdists <- matrix(NA, nrow = 10, ncol = 3000)

for (i in 1:10){
  store_climdists[i,] <- draws_s$alpha + draws_s$`beta[4]` * seq_sc_climdists[i] +
    draws_s$`beta[5]` * seq_sc_climdists[i]^2 
}

plogis(apply(store_climdists, 1, quantile, c(0.025, 0.5, 0.975))) -> sum_stat_climdists

tibble(lower = sum_stat_climdists[1,],
       upper = sum_stat_climdists[3,],
       survived = sum_stat_climdists[2,],
       clim_dist = seq_climdists) %>% 
  ggplot(aes(x = clim_dist, y = survived)) +
  geom_jitter(data = cg_model, alpha = 0.3, height = 0.05, size = 0.2, shape = 16, width = 0.0005) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.8, fill = "gray47") +
  geom_line(linewidth = 1, color = "black") +
  theme_classic(base_size = 16) +
  labs(y = "P(survival)",
       x = "climate distance") -> climdist_plot_s

## Subplot for fecundity ####
seq_sc_climdist <- seq(min(cg_model_fecun$clim_dist_sc),
                       max(cg_model_fecun$clim_dist_sc),
                       length.out = 10)
seq_climdist <- seq(min(cg_model_fecun$clim_dist),
                    max(cg_model_fecun$clim_dist),
                    length.out = 10)

store_climdist <- matrix(NA, nrow = 10, ncol = 3000)

for (i in 1:10){
  store_climdist[i,] <- draws$alpha + draws$`beta[4]` * seq_sc_climdist[i] +
    draws$`beta[5]` * seq_sc_climdist[i]^2 
}

apply(store_climdist, 1, quantile, c(0.025, 0.5, 0.975)) -> sum_stat_climdist

cg_model_fecun$ln_seed_count <- log(cg_model_fecun$seed_count)

tibble(lower = sum_stat_climdist[1,],
       upper = sum_stat_climdist[3,],
       ln_seed_count = sum_stat_climdist[2,],
       clim_dist = seq_climdist) %>% 
  ggplot(aes(x = clim_dist, y = ln_seed_count)) +
  geom_point(data = cg_model_fecun %>% filter(seed_count > 0), alpha = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.8, fill = "gray47") +
  geom_line(linewidth = 1, color = "black") +
  theme_classic(base_size = 16) +
  labs(y = "ln(seed count)",
       x = "climate distance") +
  ylim(0,10)-> climdist_plot

## Fitness subplot ####
fitness <- log(exp(sum_stat_climdist) * sum_stat_climdists)

tibble(lower = fitness[1,],
       upper = fitness[3,],
       ln_fitness = fitness[2,],
       clim_dist = seq_climdist) %>% 
  ggplot(aes(x = clim_dist, y = ln_fitness)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "gray47") +
  geom_line(linewidth = 1, color = "black") +
  theme_classic(base_size = 16) +
  labs(y = "ln(fitness)",
       x = "climate distance") -> climdist_plot_fitness

climdist_plot_s / climdist_plot / climdist_plot_fitness +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = "a") -> local_adaptation

png("figs/Fig5_localadapt.png", height = 12, width = 5.5, res = 300, units = "in")
local_adaptation
dev.off()
