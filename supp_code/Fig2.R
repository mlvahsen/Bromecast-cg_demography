# Makes plot of soil probe data (Figure 2)

# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork); library(lubridate); 
library(loo)
# Source data for modeling 
source("supp_code/data_prep.R")

# Make data set of just plants that survived and reproduced
cg_model$survived <- ifelse((cg_model$inflor_mass) > 0, 1, 0)

# Make subset of data that survived
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

# Make site year gravel index
cg_model_fecun$site_year_gravel <- paste(cg_model_fecun$site, cg_model_fecun$year, cg_model_fecun$albedo, sep = "_")

cg_model_fecun %>% 
  dplyr::select(site_year_gravel,temp_fecun, vwc_avg) %>%
  distinct() %>% 
  mutate(site = substr(site_year_gravel, 1, 2),
         year = paste(20, substr(site_year_gravel, 6, 7), sep = ""),
         site_year = paste(site, year, sep = "\n"),
         gravel = substr(site_year_gravel, 9, 13)) -> plot_data

plot_data %>% 
  ggplot(aes(x = temp_fecun, y = vwc_avg)) +
  geom_point(size = 5, aes(color = site, fill = gravel, shape = year), stroke = 2, alpha = 0.8) +
  scale_fill_manual(values = c("black", "white")) +
  scale_color_manual(values = c("#44AA99", "#AA4499", "#332288", "#6699CC")) +
  scale_shape_manual(values = c(21,22)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw(base_size = 16) +
  labs(y = "average daily soil moisture\n(volumetric water content)",
       x = "average daily temperature (°C)") -> xy_plot

png("figs/Fig2_microclimate.png", height = 5.6, width = 7.5, res = 300, units = "in")
xy_plot 
dev.off()
