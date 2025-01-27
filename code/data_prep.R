## Preliminaries ####

# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork)

# Read in full common garden data set
cg_data <- read_csv("https://raw.githubusercontent.com/pbadler/bromecast-data/main/gardens/deriveddata/cg_fullData_withFlags.csv")

## Data cleaning ####
cg_data %>%
  filter(!grepl("smut", note_standard_phen) &
           !grepl("smut", note_standard_harvest) &
           !grepl("herbivory", note_standard_phen) &
           !grepl("herbivory", note_standard_harvest) &
           !grepl("physical_damage", note_standard_phen) &
           !grepl("physical_damage", note_standard_harvest) &
           !grepl("seed_drop", note_standard_phen) &
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
plot(calib_mod) # Not horrible
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
  mutate(new_neighbors = case_when(density == "low" & possible_neighbors == 3 ~ prop_survived + neighbors,
                                   density == "low" & possible_neighbors == 2 ~ prop_survived * 2 + neighbors,
                                   density == "low" & possible_neighbors == 1 ~ prop_survived * 3 + neighbors,
                                   density == "high" & site != "WI" & possible_neighbors < 80 ~ prop_survived * (80-possible_neighbors) + neighbors,
                                   density == "high" & site == "WI" & possible_neighbors < 90 ~ prop_survived * (90-possible_neighbors) + neighbors,
                                   density == "low" & possible_neighbors > 3 ~ neighbors)) -> cg_model
##################################
## THIS IS WHERE I STOPPED ON 1/27
##################################

## Calculate climate of origin difference for each genotype × site combination ####

source("modeling/code/climate_data.R")

# Merge together climate distances with rest of data
merge(for_model, clim_dists) -> for_model

# Create unique plot ID
for_model %>% 
  mutate(plot_unique = paste(site, block, plot, sep = "_")) -> for_model



## Get weather data for gravel and site for each year ####
vwc <- read_csv("modeling/data/weather_data/dailyVWCdata_allgardens_allyears.csv") %>% clean_names()
temp <- read_csv("modeling/data/weather_data/dailytempdata_allgardens_allyears.csv") %>% clean_names()

## Calculate avg vwc for each gravel × site × year combination ####
vwc %>% 
  mutate(date = mdy(date),
         site = case_when(site == "Sheep Station" ~ "SS",
                          site == "Boise Low" ~ "WI",
                          site == "Boise High" ~ "BA",
                          site == "Cheyenne" ~ "CH"),
         gravel = ifelse(graveltrt == "Black", "black", "white")) %>% 
  dplyr::select(date, site, gravel, vwc_avg) -> vwc

vwc %>% 
  filter( (date > as_date("2022-03-10") & date < as_date("2022-05-15")) |
            date > as_date("2023-03-10") & date < as_date("2023-05-15")) %>% 
  mutate(year = year(date)) %>% 
  group_by(site, gravel, year) %>% 
  summarize(vwc_avg = mean(vwc_avg)) %>% 
  ungroup() -> vwc_sub

## Calculate avg temp for each gravel × site × year combination ####
temp %>% 
  mutate(date = mdy(date),
         site = case_when(site == "Sheep Station" ~ "SS",
                          site == "Boise Low" ~ "WI",
                          site == "Boise High" ~ "BA",
                          site == "Cheyenne" ~ "CH"),
         gravel = ifelse(graveltrt == "Black", "black", "white")) %>% 
  dplyr::select(date, site, gravel, avg_temp_c) -> temp

temp %>% 
  filter(site == "SS") %>% 
  filter(date > as_date("2023-04-24") & date < as_date("2023-08-04")) %>% 
  group_by(gravel) %>% 
  summarize(mean = mean(avg_temp_c))
  

temp %>% 
  filter( (date > as_date("2022-01-01") & date < as_date("2022-05-15")) |
            date > as_date("2023-01-01") & date < as_date("2023-05-15")) %>% 
  mutate(year = year(date),
         avg_temp_c = ifelse(avg_temp_c < 0, 0, avg_temp_c)) %>% 
  group_by(site, gravel, year) %>% 
  summarize(temp_avg_surv = mean(avg_temp_c)) %>% 
  ungroup() -> temp_surv

temp %>% 
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

# Merge together with data
merge(for_model, site_info) -> for_model

## Remove all intermediate data objects ####

rm(list=setdiff(ls(), "for_model"))

## Extra plots ####

# png("~/Documents/Research/USU/Presentations/prop_neighbors.png", height = 8.7, width = 8.7, res = 300, units = "in")
# for_model %>%
#   filter(is.na(note_standard)) %>% 
#   filter(seed_count > 0) %>%
#   ggplot(aes(x = sqrt(new_neighbors), y = log(seed_count),
#              color = albedo)) +
#   facet_grid(site ~ year) +
#   geom_point(alpha = 0.3) +
#   scale_color_manual(values = c("black", "pink")) +
#   theme_bw(base_size = 16) +
#   geom_smooth(method = "lm", aes(fill = albedo), se = F) +
#   labs(y = "ln(number of seeds)",
#        x = "sqrt(number of neighbors)") +
#   scale_fill_manual(values = c("black", "pink"))
# dev.off()
# 
# png("~/Documents/Research/USU/Presentations/pc_dist.png", height = 8.5, width = 10, res = 300, units = "in")
# for_model %>% 
#   filter(seed_count > 0) %>% 
#   ggplot(aes(x = pc_dist_score, y = log(seed_count))) +
#   geom_point() +
#   facet_wrap(~site) +
#   geom_smooth(method = "lm") +
#   theme_bw(base_size = 16) +
#   labs(x = "PC distance from common garden",
#        y = "ln(number of seeds)") -> a
# 
# for_model %>% 
#   ggplot(aes(x = pc_dist_score,
#              y = ifelse(for_model$survived == "Y", 1, 0) )) +
#   geom_point() +
#   facet_wrap(~site) +
#   geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   theme_bw(base_size = 16) +
#   labs(x = "PC distance from common garden",
#        y = "P(survival)") -> b
# 
# library(patchwork)
# a / b
# dev.off()

# colors <- c("#44AA99", "#AA4499", "#332288", "#6699CC")
# 
# site_info %>% 
#   ggplot(aes(x = year, y = vwc_avg, color = site, linetype = gravel)) +
#   geom_line(size = 2) +
#   geom_point(size = 6) +
#   theme_bw(base_size = 14) +
#   scale_x_continuous(breaks = 2022:2023) +
#   ylab("average vwc") +
#   ggtitle("soil moisture") +
#   scale_color_manual(values = colors)->vwc_plot
# 
# site_info %>% 
#   ggplot(aes(x = year, y = temp_surv, color = site, linetype = gravel)) +
#   geom_line(size = 2) +
#   geom_point(size = 6) +
#   theme_bw(base_size = 14) +
#   scale_x_continuous(breaks = 2022:2023) +
#   ylab("average temp (jan - may)") +
#   ggtitle("survival temp") +
#   ylim(2,15)+
#   scale_color_manual(values = colors)-> temp_surv_plot
# 
# site_info %>% 
#   ggplot(aes(x = year, y = temp_fecun, color = site, linetype = gravel)) +
#   geom_line(size = 2) +
#   geom_point(size = 6) +
#   theme_bw(base_size = 14) +
#   scale_x_continuous(breaks = 2022:2023) +
#   ylab("average temp (jan - may)") +
#   ggtitle("fecundity temp") +
#   scale_color_manual(values = colors)+
#   ylim(2,15)-> temp_fecun_plot
# 
# png("modeling/figs/weather_cov.png", height = 5, width = 14, res = 300, units = "in")
# temp_surv_plot + temp_fecun_plot + vwc_plot + plot_layout(guides = "collect")
# dev.off()
