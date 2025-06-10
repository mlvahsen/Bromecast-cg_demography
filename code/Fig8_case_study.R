# Plot predicted fitness for case study genotypes to see predicted effects
# (Shriver01 [Joshua Tree, CA], Chynoweth05 [Vernal, UT], AAFC09
# [Summerland, BC])

## Preliminaries ####
# Read in fecundity model output
fit <- as_cmdstan_fit(files = c("outputs/demo_model_fecun-202506091642-1-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-2-6d200d.csv",
                                "outputs/demo_model_fecun-202506091642-3-6d200d.csv"))

fit_s <- as_cmdstan_fit(files = c("outputs/demo_model_surv_noncenter-202506091738-1-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-2-51da52.csv",
                                  "outputs/demo_model_surv_noncenter-202506091738-3-51da52.csv"))

# Read in data
cg_model <- read_csv("data/cg_model_data.csv")

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
  theme(legend.position = "top") +
  scale_color_manual(values = c("#fc8d62", "#66c2a5", "#8da0cb")) -> case_1

png("figs/Fig8_casestudy.png", height = 5.4, width = 5.8, res = 300, units = "in")
case_1
dev.off()

# Get effect size for Joshua Tree (calculate median instead of mean bc very skewed)
median(exp(lnfit_gen91_high) / exp(lnfit_gen91_low))
# 3.598363 times more seeds on the original scale

# Get effect size for Joshua Tree (calculate median instead of mean bc very skewed)
median(exp(lnfit_gen5_low) / exp(lnfit_gen5_high))
# 2.791374 times more seeds on the original scale
