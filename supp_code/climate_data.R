# Load libraries
library(tidyverse); library(here); library(stringi)
library(prism); library(raster); library(geosphere); library(sf);
library(dismo); library(factoextra); library(ggfortify);
library(daymetr); library(lubridate)

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
            tmax_avg_mean = mean(tmax_avg),
            tmin_avg_mean = mean(tmin_avg)) %>% 
  ungroup() -> for_bioclim

# Calculate bioclimatic variables using for loop
store_bioclim <- matrix(NA, nrow = length(unique(for_bioclim$site)), ncol = 19)
unique_sites <- unique(for_bioclim$site)

for (i in 1:nrow(store_bioclim)){
  test <- for_bioclim %>% filter(site == unique_sites[i])
  
  store_bioclim[i,] <- biovars(prec = test$precip_total_mean,
                               tmin = test$tmin_avg_mean,
                               tmax = test$tmax_avg_mean)
}

colnames(store_bioclim) <- paste("bioclim", 1:19, sep = "_")

cbind(site_code = unique_sites, as_tibble(store_bioclim)) -> df_bioclim

# Remove common garden sites
`%notin%` <- Negate(`%in%`)

df_bioclim %>% 
  filter(site_code %notin% c("SS", "CH", "BA", "WI")) -> df_bioclim_source

# Repeat process for site years for common garden
# Write csv to current location
write_csv(gps %>% filter(site %in% c("SS", "BA", "WI", "CH")), "supp_data/cg_sites.csv")

# Get daymet data for coordinates of interest
df_batch_cg <- download_daymet_batch(file_location = "supp_data/cg_sites.csv",
                                  start = 2021,
                                  end = 2023,
                                  internal = TRUE,
                                  simplify = TRUE)

# Filter to only precip and temperature data and convert to month day
df_batch_cg %>% 
  filter(measurement %in% c("prcp..mm.day.", "tmax..deg.c.", "tmin..deg.c.")) %>% 
  mutate(date_= as.Date(yday-1, origin=paste0(year, "-01-01")), 
         month= strftime(date_, "%m"), 
         day=strftime(date_,"%d")) -> df_filt_cg

df_filt_cg %>% 
  filter(measurement == "prcp..mm.day.") %>% 
  group_by(site, year, month) %>% 
  summarize(precip_total = sum(value)) %>% 
  ungroup() -> precip_summary_cg

df_filt_cg %>% 
  filter(measurement == "tmax..deg.c.") %>% 
  group_by(site, year, month) %>% 
  summarize(tmax_avg = mean(value)) %>% 
  ungroup() -> tmax_summary_cg

df_filt_cg %>% 
  filter(measurement == "tmin..deg.c.") %>% 
  group_by(site, year, month) %>% 
  summarize(tmin_avg = mean(value)) %>% 
  ungroup() -> tmin_summary_cg

# Bring together all climate data
weathernorm <- cbind(precip_summary_cg, tmax_avg = tmax_summary_cg$tmax_avg,
                  tmin_avg = tmin_summary_cg$tmin_avg)

# For weather data, we want to create hydrological year datasets
weathernorm %>% 
  mutate(hydro_year = case_when(as.numeric(month) %in% 1:9 & year == 2022 ~ 2022,
                                as.numeric(month) %in% 10:12 & year == 2021 ~ 2022,
                                as.numeric(month) %in% 1:9 & year == 2023 ~ 2023,
                                as.numeric(month) %in% 10:12 & year == 2022 ~ 2023,
                                T ~ NA)) %>% 
  arrange(hydro_year, month) %>% 
  dplyr::select(-year, year = hydro_year) %>% 
  filter(complete.cases(year)) %>% 
  mutate(site_year = paste(site, year, sep = "_")) %>% 
  dplyr::select(-site, -year) -> for_bioclim_cg

# Calculate bioclimatic variables using for loop
store_bioclim_cg <- matrix(NA, nrow = length(unique(for_bioclim_cg$site_year)), ncol = 19)
unique_sites_cg <- unique(for_bioclim_cg$site_year)

for (i in 1:nrow(store_bioclim_cg)){
  test <- for_bioclim_cg %>% filter(site_year == unique_sites_cg[i])
  
  store_bioclim_cg[i,] <- biovars(prec = test$precip_total,
                               tmin = test$tmin_avg,
                               tmax = test$tmax_avg)
}

colnames(store_bioclim_cg) <- paste("bioclim", 1:19, sep = "_")

cbind(site_code = unique_sites_cg, as_tibble(store_bioclim_cg)) -> df_bioclim_cg


# Bring together common garden weather and genotype climate norms
rbind(df_bioclim_cg, df_bioclim_source) -> df_bioclim_all

# Center and scale all bioclimatic variables
scale(df_bioclim_all[,2:ncol(df_bioclim_all)])[,1:19] -> scaled_bioclim_all

# Save means and sds for scaling
attr(scale(df_bioclim_all[,2:ncol(df_bioclim_all)]), "scaled:center") -> means
attr(scale(df_bioclim_all[,2:ncol(df_bioclim_all)]), "scaled:scale") -> sds

write_csv(as_tibble(means), "outputs/bioclim_means.csv")
write_csv(as_tibble(sds), "outputs/bioclim_sds.csv")
# Remove all ancillary data sets
rm(list=setdiff(ls(), c("scaled_bioclim_all", "df_bioclim_all", "cg_model")))

## Run PCA ####
# Run PCA
pca_out <- prcomp(scaled_bioclim_all)

# Save PCA output for predictions later
write_rds(pca_out, "outputs/pca_out.rds")
write_rds(df_bioclim_all, "outputs/pca_df.rds")

# Get percent explained by each PC axis
round(pca_out$sdev^2 / sum(pca_out$sdev^2),3) -> perc_explained

# Figure out which variables are contributing most to PC1
sort(get_pca_var(pca_out)$contrib[,1], decreasing = T)
# Temperature and precipitation 

# Figure out which variables are contributing most to PC1
sort(get_pca_var(pca_out)$contrib[,2], decreasing = T)
# Temperature annual range, temperature seasonality, precip of coldest quarter,
# temperature of wettest quarter


## PCA of climate and weather data ####

## CODE FOR MAKING PLOT ABOUT PC2 ####
# all_clim %>% 
#   mutate(code = ifelse(genotype %in% c("SS_22", "SS_23", "BA_22", "WI_22", "WI_23",
#                                 "CH_22", "CH_23"), "CG weather", "genotype climate")) -> all_clim
# 
# png("~/Desktop/PC_axes.png", height = 5, width = 5, res = 300, units = "in")
# autoplot(pca_out, data = all_clim, color = "code", size = "code") +
#   theme_bw(base_size = 16) +
#   scale_size_manual(values = c(4,2)) +
#   scale_color_manual(values = c("orange", "black")) +
#   theme(legend.position = "top") +
#   labs(color = "",
#        size = "")
# dev.off()

# More NEGATIVE values are hotter and drier; more POSITIVE values are cooler and
# wetter

# Bind PC axis data with original data
cbind(df_bioclim_all, pca_out$x) -> bioclim_pc

# Calculate PC distance between genotypes/site × year combos
pc1 <- bioclim_pc %>% rename(pc1 = PC1) %>% pull(pc1)

# Figure out which columns are for the different sites
which(bioclim_pc$site_code %in% c("BA_2022", "CH_2022",
                                 "SS_2022", "WI_2022",
                                 "BA_2023", "CH_2023",
                                 "SS_2023", "WI_2023")) -> site_ids

# Calculate absolute distance between PC1 for source climate and PC1 for weather
`%notin%` <- Negate(`%in%`)

bioclim_pc %>%
  filter(site_code %notin% c("BA_2022", "CH_2022","SS_2022", "WI_2022",
                            "BA_2023", "CH_2023","SS_2023", "WI_2023")) %>% 
  mutate(BA_22 = PC1 - pc1[site_ids[1]],
         CH_22 = PC1 - pc1[site_ids[2]],
         SS_22 = PC1 - pc1[site_ids[3]],
         WI_22 = PC1 - pc1[site_ids[4]],
         BA_23 = PC1 - pc1[site_ids[5]],
         CH_23 = PC1 - pc1[site_ids[6]],
         SS_23 = PC1 - pc1[site_ids[7]],
         WI_23 = PC1 - pc1[site_ids[8]]) %>%
  dplyr::select(site_code, BA_22, CH_22, SS_22, WI_22,
                BA_23, CH_23, SS_23, WI_23, PC1) -> genotypes_pcs

# Get conversion information between genotype numeric code and alpha-numeric
# code
gps <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/refs/heads/main/gardens/deriveddata/BioclimateOfOrigin_AllGenotypes.csv") 

gps %>% 
  dplyr::select(site_code, genotype) -> genotype_codes

# Merge genotype codes back with calculated PC distances and PC1 and PC2
merge(genotype_codes, genotypes_pcs) -> merged_genotype_dists

merged_genotype_dists %>%
  gather(key = site_year, value = clim_dist, `BA_22`:`WI_23`) %>%
  mutate(year = 2000 + parse_number(site_year),
         site = substr(site_year, 1, 2)) %>%
  dplyr::select(site, year, genotype, clim_dist, PC1) -> clim_dists

# Remove BA 2023 because we don't have data from that year x site
clim_dists %>% 
  filter(site != "BA" | year != 2023) -> clim_dists

# Remove all data except clim_dists
rm(list=setdiff(ls(), c("clim_dists", "cg_model")))
  