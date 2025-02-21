## Preliminaries ####

# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork)

# Read in full common garden data set
cg_data <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/cg_fullData_withFlags.csv")

## Data cleaning ####
cg_data %>%
  filter(!grepl("seed_drop", note_standard_phen) &
           !grepl("seed_drop", note_standard_harvest) &
           !grepl("wrong_spp", note_standard_phen) &
           !grepl("wrong_spp", note_standard_harvest) &
           !grepl("missing", note_standard_harvest)) -> cg_clean

# Get true positives (alive at last phenology check before harvest and
# successfully harvested) and true negatives (not alive at last phenology check
# before harvest and not harvested)
cg_clean %>% 
  filter(last_phen_status == "Y" & (inflor_mass + veg_mass) > 0) -> cg_clean_pos
cg_clean %>% 
  filter(last_phen_status == "N" & (inflor_mass + veg_mass) == 0) -> cg_clean_neg
# See how many additional observations this drops out
(nrow(cg_clean_pos) + nrow(cg_clean_neg)) / nrow(cg_clean)
# Only takes out an additional 5%
rbind(cg_clean_pos, cg_clean_neg) -> cg_model
# This is a very conservative cleaning of the data

# Remove all intermediate datasets
rm(list=setdiff(ls(), "cg_model"))

# Fit a relationship between seed count and inflorescence mass for 2022 data. We
# will use this right now to calculate seed count for 2023 data.
cg_model %>% 
  filter(year == 2022 & inflor_mass > 0 & seed_count_total > 0) -> cg_calibrate

calib_mod <- lm(log(seed_count_total) ~ log(inflor_mass), data = cg_calibrate)
# plot(calib_mod) # Not horrible
calib_coefs <- as.numeric(coef(calib_mod))

# Estimate seed count for 2023 before fitting model (make this within the model
# later??)
cg_model %>% 
  mutate(seed_count = ifelse(year == 2023 & inflor_mass > 0,
                             round(exp(calib_coefs[1] + calib_coefs[2]*log(inflor_mass))),
                             seed_count_total)) -> cg_model 

# Set seed count to be 0 for any observations with NA seed count
cg_model[is.na(cg_model$seed_count), "seed_count"] <- 0

# Make a column of "survived to flowering"
cg_model %>% 
  mutate(survived_to_flower = ifelse(inflor_mass > 0, "Y", "N")) -> cg_model

## Calculate neighborhood density for each plant ####

# Set possible number of neighbors for each location in high density
cg_model %>%
  mutate(plot_unique = paste(site, block, plot, sep = "_")) -> cg_model

cg_model$possible_neighbors <- NULL
cg_model$neighbors <- NULL
cg_model$prop_neighbors <- NULL

for(i in 1:nrow(cg_model)){
  
  if(cg_model$density[i] == "lo"){
    cg_model[i,] %>% 
      dplyr::select(x, y) %>% 
      mutate(x_new = x + 1,
             x_new2 = x - 1,
             y_new = y + 1,
             y_new2 = y - 1) -> search_coords
    
    cg_model %>% 
      filter(plot_unique == cg_model$plot_unique[i]) %>% 
      filter(x == search_coords$x_new & y == search_coords$y |
               x == search_coords$x_new2 & y == search_coords$y |
               x == search_coords$x & y == search_coords$y_new  |
               x == search_coords$x & y == search_coords$y_new2 ) -> possible_neighbors
  }else{
      expand.grid(x = cg_model[i,]$x + -5:5, y = cg_model[i,]$y + -5:5) -> search_coords
      
    # Filter out search coords that are not within circle using distance matrix
    distances <- as.matrix(dist(cbind(search_coords$x, search_coords$y)))
    focal_coords <- which(search_coords$x == cg_model$x[i] & search_coords$y == cg_model$y[i])
    search_coords <- search_coords %>% 
      mutate(dist = distances[focal_coords,]) %>% 
      filter(dist <= 5)
    
    cg_model %>% 
        filter(plot_unique == cg_model$plot_unique[i]) %>% 
        filter(x %in% search_coords$x & y %in% search_coords$y) %>% 
        filter(x != cg_model$x[i] | y != cg_model$y[i]) -> possible_neighbors
  }

  cg_model[i, "possible_neighbors"] <- nrow(possible_neighbors)
  cg_model[i, "neighbors"] <- nrow(possible_neighbors %>% filter(survived_to_flower == "Y"))

}

## Adjust for edge effects ####

# Get proportion that survived for each plot
cg_model %>% 
  mutate(w = ifelse(survived_to_flower == "Y", 1, 0)) %>% 
  group_by(plot_unique) %>% 
  summarize(prop_survived = sum(w)/n()) %>% 
  ungroup() -> plot_survival

merge(cg_model, plot_survival) -> cg_model

cg_model %>% 
  mutate(new_neighbors = case_when(density == "lo" & possible_neighbors == 3 ~ prop_survived + neighbors,
                                   density == "lo" & possible_neighbors == 2 ~ prop_survived * 2 + neighbors,
                                   density == "lo" & possible_neighbors == 1 ~ prop_survived * 3 + neighbors,
                                   # for 2023 there were less possible neighbors because there were less plants (WI had up to 90, all other sites up to 80)
                                   density == "hi" & site != "WI" & possible_neighbors < 80 ~ prop_survived * (80-possible_neighbors) + neighbors,
                                   density == "hi" & site == "WI" & possible_neighbors < 90 ~ prop_survived * (90-possible_neighbors) + neighbors,
                                   density == "lo" & possible_neighbors > 3 ~ neighbors)) -> cg_model
##################################
## THIS IS WHERE I STOPPED ON 1/27
##################################

## Calculate climate of origin difference for each genotype × site combination ####
source("code/climate_data.R")

# Merge together climate distances with rest of data (right now this drops out
# all Canadian genotypes)
merge(cg_model, clim_dists) -> cg_model

# Create unique plot ID
cg_model %>% 
  mutate(plot_unique = paste(site, block, plot, sep = "_")) -> cg_model

## Get weather data for gravel and site for each year ####
vwc <- read_csv("supp_data/weather_data/dailyVWCdata_allgardens_allyears.csv") %>% clean_names()
temp <- read_csv("supp_data/weather_data/dailytempdata_allgardens_allyears.csv") %>% clean_names()

## Calculate avg vwc for each gravel × site × year combination ####

vwc$date <- mdy(vwc$date)

# Use Mar - May window for consistent data continuity across sites
vwc %>% 
  filter( (date > as_date("2022-03-10") & date < as_date("2022-05-15")) |
            date > as_date("2023-03-10") & date < as_date("2023-05-15")) %>% 
  mutate(year = year(date)) %>%
  rename(gravel = graveltrt) %>% 
  group_by(site, gravel, year) %>% 
  summarize(vwc_avg = mean(vwc_avg)) %>% 
  ungroup() -> vwc_sub

## Calculate avg temp for each gravel × site × year combination ####

# Use Jan - May window for survival part of the model
temp %>% 
  mutate(date = mdy(date)) %>% 
  rename(gravel = graveltrt) %>% 
  filter( (date > as_date("2022-01-01") & date < as_date("2022-05-15")) |
            date > as_date("2023-01-01") & date < as_date("2023-05-15")) %>% 
  mutate(year = year(date),
         avg_temp_c = ifelse(avg_temp_c < 0, 0, avg_temp_c)) %>% 
  group_by(site, gravel, year) %>% 
  summarize(temp_avg_surv = mean(avg_temp_c)) %>% 
  ungroup() -> temp_surv

# Use Mar - May window for fecundity part of the model
temp %>% 
  mutate(date = mdy(date)) %>% 
  rename(gravel = graveltrt) %>% 
  filter( (date > as_date("2022-03-01") & date < as_date("2022-05-15")) |
            date > as_date("2023-03-01") & date < as_date("2023-05-15")) %>% 
  mutate(year = year(date),
         avg_temp_c = ifelse(avg_temp_c < 0, 0, avg_temp_c)) %>% 
  group_by(site, gravel, year) %>% 
  summarize(temp_avg_fecun = mean(avg_temp_c)) %>% 
  ungroup() -> temp_fecun

# Bring all data climate data together
site_info <- cbind(vwc_sub,
                   temp_surv = temp_surv$temp_avg_surv,
                   temp_fecun = temp_fecun$temp_avg_fecun) %>% 
  rename(albedo = gravel)

# Fix gravel and site variables
site_info %>% 
  mutate(site = case_when(site == "Sheep Station" ~ "SS",
                          site == "Boise High" ~ "BA",
                          site == "Boise Low" ~ "WI",
                          site == "Cheyenne" ~ "CH"),
         albedo = ifelse(albedo == "Black", "black", "white")) -> site_info

# Merge together with data
merge(cg_model, site_info) -> cg_model

## Remove all intermediate data objects ####
rm(list=setdiff(ls(), "for_model"))
