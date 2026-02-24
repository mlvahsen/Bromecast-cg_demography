# Load libraries
library(tidyverse); library(here); library(stringi)
library(prism); library(raster); library(geosphere); library(sf);
library(dismo); library(factoextra); library(ggfortify);
library(daymetr); library(lubridate); library(scales)

## Calculate prediction matrices ####

# Read GPS data for all genotypes from main Bromecast repository
gps <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/refs/heads/main/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv") %>%
  dplyr::select(site = site_code, lat, lon)

# Write csv to current location
write_csv(gps, "supp_data/gps_sites.csv")

# Get daymet data for coordinates of interest
df_batch <- download_daymet_batch(file_location = "supp_data/gps_sites.csv",
                                  start = 1991,
                                  end = 2020,
                                  internal = TRUE,
                                  simplify = TRUE)

# Filter to only precip and temperature data and convert to month day
df_batch %>%
  filter(measurement %in% c("prcp..mm.day.", "tmax..deg.c.", "tmin..deg.c.")) %>%
  mutate(date_= as.Date(yday-1, origin=paste0(year, "-01-01")),
         month= strftime(date_, "%m"),
         day=strftime(date_,"%d")) -> df_filt

df_filt %>%
  filter(measurement == "prcp..mm.day.") %>%
  group_by(site, year, month) %>%
  summarize(precip_total = sum(value)) %>%
  ungroup() -> precip_summary

df_filt %>%
  filter(measurement == "tmax..deg.c.") %>%
  group_by(site, year, month) %>%
  summarize(tmax_avg = mean(value)) %>%
  ungroup() -> tmax_summary

df_filt %>%
  filter(measurement == "tmin..deg.c.") %>%
  group_by(site, year, month) %>%
  summarize(tmin_avg = mean(value)) %>%
  ungroup() -> tmin_summary

# Bring together all climate data
climnorm <- cbind(precip_summary, tmax_avg = tmax_summary$tmax_avg,
                  tmin_avg = tmin_summary$tmin_avg)

# For climate norm, we want the average for each month across all of the years
climnorm %>%
  group_by(site, month) %>%
  summarize(precip_total_mean = mean(precip_total),
            # Add 2 degrees here
            tmax_avg_mean = mean(tmax_avg + 2),
            tmin_avg_mean = mean(tmin_avg + 2)) %>%
  ungroup() -> for_bioclim

# Calculate bioclimatic variables using for loop for REAL SITES + 2 degrees Celsius
store_bioclim_future <- matrix(NA, nrow = length(unique(for_bioclim$site)), ncol = 19)
unique_sites <- unique(for_bioclim$site)

for (i in 1:nrow(store_bioclim_future)){
  test <- for_bioclim %>% filter(site == unique_sites[i])
  
  store_bioclim_future[i,] <- biovars(prec = test$precip_total_mean,
                                      tmin = test$tmin_avg_mean,
                                      tmax = test$tmax_avg_mean)
}

colnames(store_bioclim_future) <- paste("bioclim", 1:19, sep = "_")

cbind(site_code = unique_sites, as_tibble(store_bioclim_future)) -> df_bioclim_future

`%notin%` <- Negate(`%in%`)

df_bioclim_future %>%
  filter(site_code %notin% c("BA_2022", "CH_2022",
                             "SS_2022", "WI_2022",
                             "BA_2023", "CH_2023",
                             "SS_2023", "WI_2023", "SS", "CH", "WI", "BA")) -> df_bioclim_source_future

# Rename sites so we can tell which ones are from the future
df_bioclim_source_future %>%
  mutate(site_code = paste(site_code, "future", sep = "_")) -> df_bioclim_source_future

# Center and scale all bioclimatic variables
means <- read_csv("outputs/bioclim_means.csv") %>% t()
sds <- read_csv("outputs/bioclim_sds.csv") %>% t()

sweep(df_bioclim_source_future[,2:20], 2, means, FUN = "-") -> mat_centered
sweep(mat_centered, 2, sds, FUN = "/") -> mat_scaled

# Run PCA
# Load fitted pca
pca_out <- read_rds("outputs/pca_out.rds")
pca_site_codes <- read_rds("outputs/pca_df.rds")

predict(pca_out, mat_scaled)[,1] -> PC1_future

# Get future PC1 value for each site
tibble(site_code = df_bioclim_source_future$site_code,
       PC1_future = PC1_future) %>%
  mutate(site_code = sub("_future$", "", site_code)) -> future_pc1

# Get original PC1 value for each site
tibble(site_code = pca_site_codes$site_code,
       PC1_source = pca_out$x[,1]) %>%
  merge(future_pc1) %>%
  # Get climate difference by subtracting future climate (source climate + 2) and
  # source climate
  mutate(climate_diff = PC1_source - PC1_future) -> pred_matrix
# Ok so this is the same value for each one
# 1.109024
# Get annual temperature for each site-year combo in the garden
pca_df <- read_rds("outputs/pca_df.rds")

pca_df %>%
  filter(site_code %in% c("BA_2022", "WI_2022", "SS_2022", "CH_2022",
                          "WI_2023", "SS_2023", "CH_2023")) %>%
  dplyr::select(site_code, bioclim_1) %>%
  mutate(site = substr(site_code, 1, 2),
         year = parse_number(site_code)) -> annual_temp

cg_model %>%
  group_by(site, year) %>%
  summarize(soil_temp_fecun = mean(temp_fecun),
            soil_temp_surv = mean(temp_surv)) %>%
  ungroup() %>%
  merge(annual_temp) -> anntemp_soiltemp

mod_temp_fecun <- lm(soil_temp_fecun ~ bioclim_1, data = anntemp_soiltemp)
summary(mod_temp_fecun)# 0.8239
mod_temp_surv <- lm(soil_temp_surv ~ bioclim_1, data = anntemp_soiltemp)
summary(mod_temp_surv) # 0.8148

# Get predicted soil temp for each source site based on 30-year climate norm
pca_df %>%
  mutate(soil_temp_fecun = coef(mod_temp_fecun)[1] + coef(mod_temp_fecun)[2]*bioclim_1,
         soil_temp_surv = coef(mod_temp_surv)[1] + coef(mod_temp_surv)[2]*bioclim_1) %>%
  dplyr::select(site_code, soil_temp_fecun, soil_temp_surv) %>%
  merge(pred_matrix) %>%
  mutate(soil_temp_fecun_pred = soil_temp_fecun + 2,
         soil_temp_surv_pred = soil_temp_surv + 2) %>%
  dplyr::select(site_code, soil_temp_fecun, soil_temp_fecun_pred,
                soil_temp_surv, soil_temp_surv_pred, climate_diff) -> pred_matrix

write_csv(pred_matrix, "outputs/pred_future_matrix.csv")

cg_model <- read_csv("data/cg_model_data.csv")# Updated 3 April 2025

# Make data set of just plants that survived and reproduced
cg_model$survived <- ifelse(cg_model$seed_count > 0, 1, 0)

# Make subset of data that survived
cg_model %>%
  filter(survived == 1) -> cg_model_fecun

# Center and scale continuous variables, make factors, make survival response
# binary (from previous model fits)
clim_dist_sc_mean_f <- attr(scale(cg_model_fecun$clim_dist), "scaled:center")
clim_dist_sc_sd_f <- attr(scale(cg_model_fecun$clim_dist), "scaled:scale")

clim_dist_sc_mean_s <- attr(scale(cg_model$clim_dist), "scaled:center")
clim_dist_sc_sd_s <- attr(scale(cg_model$clim_dist), "scaled:scale")

cg_model_fecun_sc_mean <- attr(scale(cg_model_fecun$temp_fecun), "scaled:center")
cg_model_fecun_sc_sd <- attr(scale(cg_model_fecun$temp_fecun), "scaled:scale")

cg_model_surv_sc_mean <- attr(scale(cg_model$temp_surv), "scaled:center")
cg_model_surv_sc_sd <- attr(scale(cg_model$temp_surv), "scaled:scale")

genotype <- as.factor(cg_model_fecun$genotype)
genotype_id <- as.numeric(as.factor(cg_model_fecun$genotype))

# Read in GPS values for source sites
gps <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/refs/heads/main/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv")

gps %>%
  mutate(genotype = as.factor(genotype)) %>%
  dplyr::select(site_code, genotype) ->site2genotype

tibble(genotype, genotype_id) %>%
  distinct() %>%
  merge(site2genotype) %>%
  merge(pred_matrix) %>%
  arrange(genotype_id) %>%
  mutate(clim_dist_sc_f = (climate_diff - clim_dist_sc_mean_f)/clim_dist_sc_sd_f,
         clim_dist_sc2_f = clim_dist_sc_f^2,
         clim_dist_sc_s = (climate_diff - clim_dist_sc_mean_s)/clim_dist_sc_sd_s,
         clim_dist_sc2_s = clim_dist_sc_s^2,
         temp_fecun_sc = (soil_temp_fecun_pred - cg_model_fecun_sc_mean)/cg_model_fecun_sc_sd,
         temp_surv_sc = (soil_temp_surv_pred - cg_model_surv_sc_mean)/cg_model_surv_sc_sd) %>%
  dplyr::select(genotype, genotype_id, soil_temp_fecun_pred,soil_temp_surv_pred,
                climate_diff, clim_dist_sc_f, clim_dist_sc2_f,clim_dist_sc_s, clim_dist_sc2_s,
                temp_fecun_sc, temp_surv_sc) -> pred_matrix_future

pred_matrix_source <- tibble(genotype, genotype_id) %>%
  distinct() %>%
  merge(site2genotype) %>%
  merge(pred_matrix) %>%
  arrange(genotype_id) %>%
  mutate(temp_fecun_sc = (soil_temp_fecun - cg_model_fecun_sc_mean)/cg_model_fecun_sc_sd,
         temp_surv_sc = (soil_temp_surv - cg_model_surv_sc_mean)/cg_model_surv_sc_sd) %>%
  dplyr::select(genotype, genotype_id, soil_temp_fecun, temp_fecun_sc,
                soil_temp_surv, temp_surv_sc)

# Write csvs to save these for prediction
write_csv(pred_matrix_source, "outputs/predict_source_mm.csv")
write_csv(pred_matrix_future, "outputs/predict_future_mm.csv")
