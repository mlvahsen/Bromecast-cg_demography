# Plot predicted fitness for case study genotypes to see predicted effects
# (Shriver01 [Joshua Tree, CA], Chynoweth05 [Vernal, UT], AAFC09
# [Summerland, BC])

## Preliminaries ####
library(tidyverse)
library(cmdstanr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)

# Read in fecundity model output
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

## Extract covariates and model parameters ####
cg_model %>% 
  # Get scaled temp, moisture, and climate distance and squared value too
  mutate(temp_surv_sc = scale(temp_surv)[,1],
         temp_fecun_sc = scale(temp_fecun)[,1],
         vwc_avg_sc = scale(vwc_avg)[,1],
         clim_dist_sc = scale(clim_dist)[,1],
         clim_dist_sc2 = clim_dist_sc^2) %>% 
  filter(genotype %in% c(6, 98, 101)) %>% 
  filter(site == "CH" & year == 2023 & albedo == "white" | (site == "WI" & year == 2022 & albedo == "black")) %>% 
  dplyr::select(genotype, site, year, albedo, temp_surv_sc, temp_fecun_sc,
                vwc_avg_sc, clim_dist_sc, clim_dist_sc2) %>% 
  distinct() -> x_variables

cg_model_fecun %>% 
  # Get scaled temp, moisture, and climate distance and squared value too
  mutate(temp_surv_sc = scale(temp_surv)[,1],
         temp_fecun_sc = scale(temp_fecun)[,1],
         vwc_avg_sc = scale(vwc_avg)[,1],
         clim_dist_sc = scale(clim_dist)[,1],
         clim_dist_sc2 = clim_dist_sc^2) %>% 
  filter(genotype %in% c(6, 98, 101)) %>% 
  filter(site == "CH" & year == 2023 & albedo == "white" |
           (site == "WI" & year == 2022 & albedo == "black") |
           # Need to get more info for clim dist
           genotype == 101 & site == "CH" & year == 2023 & albedo == "black") %>% 
  dplyr::select(genotype, site, year, albedo, temp_surv_sc, temp_fecun_sc,
                vwc_avg_sc, clim_dist_sc, clim_dist_sc2) %>% 
  distinct() -> x_variables_f

cg_model %>% 
  # Get scaled temp, moisture, and climate distance and squared value too
  mutate(temp_surv_sc = scale(temp_surv)[,1],
         temp_fecun_sc = scale(temp_fecun)[,1],
         vwc_avg_sc = scale(vwc_avg)[,1],
         clim_dist_sc = scale(clim_dist)[,1],
         clim_dist_sc2 = clim_dist_sc^2) %>% 
  filter(genotype %in% c(6, 98, 101)) %>% 
  filter(site == "CH" & year == 2023 & albedo == "white" |
           (site == "WI" & year == 2022 & albedo == "black")) %>% 
  dplyr::select(genotype, site, year, albedo, temp_surv_sc, temp_fecun_sc,
                vwc_avg_sc, clim_dist_sc, clim_dist_sc2) %>% 
  distinct() -> x_variables_s

# Extract model draws
draws_f <- fit$draws(format = "df")
draws_s <- fit_s$draws(format = "df")

# Write function to extract coefficients
extract_coef_matrix <- function(draws, name) {
  draws %>% as_tibble() %>%
    select(starts_with(name)) %>%
    as.matrix()
}

coef_list <- list(
  alpha_f = draws_f$alpha,
  alpha_s = draws_s$alpha,
  intercepts_f = extract_coef_matrix(draws_f, "genotype_intercepts"),
  intercepts_s = extract_coef_matrix(draws_s, "genotype_intercepts"),
  beta_f = list(
    temp = draws_f$`beta[2]`,
    vwc = draws_f$`beta[3]`,
    clim = draws_f$`beta[4]`,
    clim2 = draws_f$`beta[5]`
  ),
  beta_s = list(
    temp = draws_s$`beta[2]`,
    vwc = draws_s$`beta[3]`,
    clim = draws_s$`beta[4]`,
    clim2 = draws_s$`beta[5]`
  ),
  slopes_f = list(
    temp = extract_coef_matrix(draws_f, "genotype_slopes") %>% .[, grepl(",2]", colnames(.))],
    vwc = extract_coef_matrix(draws_f, "genotype_slopes") %>% .[, grepl(",3]", colnames(.))]
  ),
  slopes_s = list(
    temp = extract_coef_matrix(draws_s, "genotype_slopes") %>% .[, grepl(",2]", colnames(.))],
    vwc = extract_coef_matrix(draws_s, "genotype_slopes") %>% .[, grepl(",3]", colnames(.))]
  )
)

## Fecundity cool and wet site ####
f_data_91 <- x_variables_f %>% filter(genotype == 101)
f_data_5 <- x_variables_f %>% filter(genotype == 6)
f_data_88 <- x_variables_f %>% filter(genotype == 98)

# Need to get scaled temperature from cg_model dataset
(cg_model %>% filter(site == "CH" & year == 2023 & albedo == "white") %>% 
    pull(temp_fecun) %>% unique() - 9.237818) / 3.223435 # -1.249972
# Need to get scaled vwc from cg_model dataset
(cg_model %>% filter(site == "CH" & year == 2023 & albedo == "white") %>% 
    pull(vwc_avg) %>% unique() - 0.1548502) / 0.03055884 # 2.048437

fecun_gen91_low <- (coef_list$alpha_f + coef_list$intercepts_f[,91]) +
  (coef_list$beta_f$temp + coef_list$slopes_f$temp[,91]) * -1.249972 +
  (coef_list$beta_f$vwc + coef_list$slopes_f$vwc[,91]) * 2.048437 +
  coef_list$beta_f$clim * f_data_91$clim_dist_sc[1] +
  coef_list$beta_f$clim2 * f_data_91$clim_dist_sc2[1]

fecun_gen5_low <- (coef_list$alpha_f + coef_list$intercepts_f[,5]) +
  (coef_list$beta_f$temp + coef_list$slopes_f$temp[,5]) * -1.249972 +
  (coef_list$beta_f$vwc + coef_list$slopes_f$vwc[,5]) * 2.048437 +
  coef_list$beta_f$clim * f_data_5$clim_dist_sc[1] +
  coef_list$beta_f$clim2 * f_data_5$clim_dist_sc2[1]

fecun_gen88_low <- (coef_list$alpha_f + coef_list$intercepts_f[,88]) +
  (coef_list$beta_f$temp + coef_list$slopes_f$temp[,88]) * -1.249972 +
  (coef_list$beta_f$vwc + coef_list$slopes_f$vwc[,88]) * 2.048437 +
  coef_list$beta_f$clim * f_data_88$clim_dist_sc[1] +
  coef_list$beta_f$clim2 * f_data_88$clim_dist_sc2[1]

## Fecundity hot & dry site ####
fecun_gen91_high <- (coef_list$alpha_f + coef_list$intercepts_f[,91]) +
  (coef_list$beta_f$temp + coef_list$slopes_f$temp[,91]) * f_data_91$temp_fecun_sc[2] +
  (coef_list$beta_f$vwc + coef_list$slopes_f$vwc[,91]) * f_data_91$vwc_avg_sc[2] +
  coef_list$beta_f$clim * f_data_91$clim_dist_sc[2] +
  coef_list$beta_f$clim2 * f_data_91$clim_dist_sc2[2]

fecun_gen5_high <- (coef_list$alpha_f + coef_list$intercepts_f[,5]) +
  (coef_list$beta_f$temp + coef_list$slopes_f$temp[,5]) * f_data_5$temp_fecun_sc[2] +
  (coef_list$beta_f$vwc + coef_list$slopes_f$vwc[,5]) * f_data_5$vwc_avg_sc[2] +
  coef_list$beta_f$clim * f_data_5$clim_dist_sc[2] +
  coef_list$beta_f$clim2 * f_data_5$clim_dist_sc2[2]

fecun_gen88_high <- (coef_list$alpha_f + coef_list$intercepts_f[,88]) +
  (coef_list$beta_f$temp + coef_list$slopes_f$temp[,88]) * f_data_88$temp_fecun_sc[2] +
  (coef_list$beta_f$vwc + coef_list$slopes_f$vwc[,88]) * f_data_88$vwc_avg_sc[2] +
  coef_list$beta_f$clim * f_data_88$clim_dist_sc[2] +
  coef_list$beta_f$clim2 * f_data_88$clim_dist_sc2[2]

## Survival cool & wet site ####
s_data_91 <- x_variables_s %>% filter(genotype == 101)
s_data_5 <- x_variables_s %>% filter(genotype == 6)
s_data_88 <- x_variables_s %>% filter(genotype == 98)

surv_gen91_low <- (coef_list$alpha_s + coef_list$intercepts_s[,91]) +
  (coef_list$beta_s$temp + coef_list$slopes_s$temp[,91]) * s_data_91$temp_surv_sc[1] +
  (coef_list$beta_s$vwc + coef_list$slopes_s$vwc[,91]) * s_data_91$vwc_avg_sc[1] +
  coef_list$beta_s$clim * s_data_91$clim_dist_sc[1] +
  coef_list$beta_s$clim2 * s_data_91$clim_dist_sc2[1]

surv_gen5_low <- (coef_list$alpha_s + coef_list$intercepts_s[,5]) +
  (coef_list$beta_s$temp + coef_list$slopes_s$temp[,5]) * s_data_5$temp_surv_sc[1] +
  (coef_list$beta_s$vwc + coef_list$slopes_s$vwc[,5]) * s_data_5$vwc_avg_sc[1] +
  coef_list$beta_s$clim * s_data_5$clim_dist_sc[1] +
  coef_list$beta_s$clim2 * s_data_5$clim_dist_sc2[1]

surv_gen88_low <- (coef_list$alpha_s + coef_list$intercepts_s[,88]) +
  (coef_list$beta_s$temp + coef_list$slopes_s$temp[,88]) * s_data_88$temp_surv_sc[1] +
  (coef_list$beta_s$vwc + coef_list$slopes_s$vwc[,88]) * s_data_88$vwc_avg_sc[1] +
  coef_list$beta_s$clim * s_data_88$clim_dist_sc[1] +
  coef_list$beta_s$clim2 * s_data_88$clim_dist_sc2[1]

## Survival hot & dry site ####

surv_gen91_high <- (coef_list$alpha_s + coef_list$intercepts_s[,91]) +
  (coef_list$beta_s$temp + coef_list$slopes_s$temp[,91]) * s_data_91$temp_surv_sc[2] +
  (coef_list$beta_s$vwc + coef_list$slopes_s$vwc[,91]) * s_data_91$vwc_avg_sc[2] +
  coef_list$beta_s$clim * s_data_91$clim_dist_sc[2] +
  coef_list$beta_s$clim2 * s_data_91$clim_dist_sc2[2]

surv_gen88_high <- (coef_list$alpha_s + coef_list$intercepts_s[,88]) +
  (coef_list$beta_s$temp + coef_list$slopes_s$temp[,88]) * s_data_88$temp_surv_sc[2] +
  (coef_list$beta_s$vwc + coef_list$slopes_s$vwc[,88]) * s_data_88$vwc_avg_sc[2] +
  coef_list$beta_s$clim * s_data_88$clim_dist_sc[2] +
  coef_list$beta_s$clim2 * s_data_88$clim_dist_sc2[2]

surv_gen5_high <- (coef_list$alpha_s + coef_list$intercepts_s[,5]) +
  (coef_list$beta_s$temp + coef_list$slopes_s$temp[,5]) * s_data_5$temp_surv_sc[2] +
  (coef_list$beta_s$vwc + coef_list$slopes_s$vwc[,5]) * s_data_5$vwc_avg_sc[2] +
  coef_list$beta_s$clim * s_data_5$clim_dist_sc[2] +
  coef_list$beta_s$clim2 * s_data_5$clim_dist_sc2[2]

## Calculate log-fitness ####

lnfit_gen91_low <- log(exp(fecun_gen91_low) * plogis(surv_gen91_low))
lnfit_gen5_low <- log(exp(fecun_gen5_low) * plogis(surv_gen5_low))
lnfit_gen88_low <- log(exp(fecun_gen88_low) * plogis(surv_gen88_low))

lnfit_gen91_high <- log(exp(fecun_gen91_high) * plogis(surv_gen91_high))
lnfit_gen5_high <- log(exp(fecun_gen5_high) * plogis(surv_gen5_high))
lnfit_gen88_high <- log(exp(fecun_gen88_high) * plogis(surv_gen88_high))

## Create plot ####
tibble(pred = c(lnfit_gen91_low, lnfit_gen5_low, lnfit_gen88_low,
                lnfit_gen91_high, lnfit_gen5_high, lnfit_gen88_high),
       genotype = rep(c("Joshua Tree, CA","Summerland, BC", "Vernal, UT",
                        "Joshua Tree, CA","Summerland, BC", "Vernal, UT"), each = 6000),
       temp = rep(c("cool & wet", "hot & dry"), each = 18000)) %>% 
  group_by(genotype, temp) %>% 
  summarize(mean = mean(pred),
            lower = quantile(pred, 0.025),
            upper = quantile(pred, 0.975)) %>% 
  ggplot(aes(x = temp, y = mean, group = genotype)) +
  geom_point(size = 4, aes(color = genotype), position = position_dodge(width = 0.3)) +
  geom_linerange(aes(x = temp, ymin = lower, ymax = upper, color = genotype),
               linewidth = 1.2, position = position_dodge2(width = 0.3)) +
  theme_classic(base_size = 16) +
  geom_line(linetype = "dashed", 
            position = position_dodge(width = 0.3),
            aes(color = genotype)) +
  labs(x = "current environment",
       y = "predicted ln(fitness)",
       color = "") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#DC267F", "#4B5099", "#785EF0")) -> case_1

# Sample GPS coordinates in Wyoming and Idaho
gps_data <- data.frame(
  site = c("Joshua Tree, CA", "Summerland, BC", "Vernal, UT"),
  latitude = c(34.06, 49.62, 40.48),
  longitude = c(-116.22, -119.8, -109.44)
)

# Get map data for Canada and the US with provinces/states
canada <- ne_states(country = "Canada", returnclass = "sf")
usa <- ne_states(country = "United States of America", returnclass = "sf")

# Combine the US states and BC for plotting
map_data_combined <- rbind(usa, canada)

ggplot() +
  geom_sf(data = map_data_combined, fill = "lightgray", color = "white",
          linewidth = 0.8) +
  coord_sf(xlim = c(-125, -107.5), ylim = c(32, 50.5), expand = FALSE) +
  geom_point(data = gps_data, aes(x = longitude, y = latitude, color = site), size = 6,
             shape = 16, stroke = 2) +
  theme_minimal(base_size = 16) +
  labs(x = "longitude", y = "latitude") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#DC267F", "#4B5099", "#785EF0")) +
  annotate(geom = "text", x = -113.22, y = 36.5, label = "Joshua Tree, CA",
           color = "#DC267F", size = 7) +
  annotate(geom = "text", x = -113.22, y = 35.6, label = "(1233 m)",
           color = "#DC267F", size = 6) +
  annotate(geom = "segment", x = -116.22, y = 34.06, xend = -115.9, yend = 36, color = "#DC267F") +
  annotate(geom = "text", x = -114, y = 41.3, label = "Vernal, UT",
           color = "#785EF0", size = 7) +
  annotate(geom = "text", x = -114, y = 40.4, label = "(1624 m)",
           color = "#785EF0", size = 6) +
  annotate(geom = "segment", x = -109.44, y = 40.48, xend = -111.5, yend = 40.7, color = "#785EF0") +
  annotate(geom = "text", x = -115.8, y = 47.62, label = "Summerland, BC",
           color = "#4B5099", size = 7) +
  annotate(geom = "text", x = -115.8, y = 46.72, label = "(920 m)",
           color = "#4B5099", size = 6) +
  annotate(geom = "segment", x = -119.8, y = 49.62, xend = -117, yend = 47.9, color = "#4B5099") -> map

png("figs/Fig6_casestudy.png", height = 5.82, width = 9.64, res = 300, units = "in")
map + case_1 + plot_annotation(tag_levels = "a")
dev.off()

# Get effect size for Joshua Tree (calculate median instead of mean bc very skewed)
median(exp(lnfit_gen91_high) / exp(lnfit_gen91_low))
# 3.598363 times more seeds on the original scale

# Get effect size for Joshua Tree (calculate median instead of mean bc very skewed)
median(exp(lnfit_gen5_low) / exp(lnfit_gen5_high))
# 2.791374 times more seeds on the original scale
